// Pairing-Based Calculator.

// TODO: Garbage collection.
#include <unistd.h>  // For getopt.

#include "pbc.h"
#include "pbc_fp.h"
#include "pbc_multiz.h"

#include "misc/darray.h"
#include "misc/symtab.h"

#include "pbc_tree.h"

#include "lex.yy.h"
#include "parser.tab.h"

int option_easy = 0;
const char *option_prompt;

char *pbc_getline(const char *prompt);

void yyerror(char *s) {
  fprintf(stderr, "%s\n", s);
}

int yyparse(void);

// Symbol table holding built-in functions and variables.
static symtab_t reserved;
// Symbol table holding user-defined variable and function names.
static symtab_t tab;

static field_t Z;
static pairing_t pairing;

// Syntax tree node.
struct tree_s {
  // Evaluates this node.
  val_ptr (*fun)(tree_ptr);
  void *data;
  // Child nodes; NULL if no children.
  darray_ptr child;
};

// Evaluates syntax tree node.
static val_ptr tree_eval(tree_ptr t) {
  return t->fun(t);
}

// The interface of a val_ptr shared amongst many val_ptr objects.
// Analog of C++ class.
struct val_type_s {
  char *name;
  void (*out_str)(FILE *, val_ptr);
  val_ptr (*eval)(val_ptr);
  val_ptr (*funcall)(val_ptr, tree_ptr);
};

// When interpreting, each node of the syntax tree recursively evaluates
// its children then returns a val_ptr.
struct val_s {
  struct val_type_s *type;
  void *data;
};

static val_ptr val_new_element(element_ptr e);
static val_ptr val_new_field(field_ptr e);

static void v_elem_out(FILE* stream, val_ptr v) {
  element_out_str(stream, 0, v->data);
}

static val_ptr v_elem_eval(val_ptr v) {
  element_ptr e = pbc_malloc(sizeof(*e));
  element_init_same_as(e, v->data);
  element_set(e, v->data);
  return val_new_element(e);
}

static void v_fun_out(FILE* stream, val_ptr v) {
  UNUSED_VAR(v);
  fprintf(stream, "function");
}

static val_ptr v_fun_call(val_ptr v, tree_ptr t) {
  val_ptr(*fun)(tree_ptr) = v->data;
  return fun(t);
}

static val_ptr v_field_cast(val_ptr v, tree_ptr t) {
  // TODO: Check args, x is an element.
  val_ptr x = tree_eval(darray_at(t->child, 0));
  element_ptr e = x->data;
  if (v->data == Z || e->field == Z) {
    // Map to/from integer.
    if (v->data == e->field) return x;
    mpz_t z;
    mpz_init(z);
    element_to_mpz(z, e);
    element_clear(e);
    element_init(e, v->data);
    element_set_mpz(e, z);
    mpz_clear(z);
  }
  return x;
}

static void v_field_out(FILE* stream, val_ptr v) {
  field_out_info(stream, v->data);
}

static val_ptr v_self(val_ptr v) {
  return v;
}

static void v_err_out(FILE* stream, val_ptr v) {
  fprintf(stream, "%s", (const char *) v->data);
}

static val_ptr v_err_call(val_ptr v, tree_ptr t) {
  UNUSED_VAR(t);
  return v;
}

static struct val_type_s
  v_elem[1]  = {{  "element",  v_elem_out, v_elem_eval, NULL }},
  v_field[1] = {{    "field", v_field_out,      v_self, v_field_cast }},
  v_fun[1]   = {{ "function",   v_fun_out,      v_self, v_fun_call }},
  v_error[1] = {{    "error",   v_err_out,      v_self, v_err_call }};

static val_ptr val_new_element(element_ptr e) {
  val_ptr v = pbc_malloc(sizeof(*v));
  v->type = v_elem;
  v->data = e;
  return v;
}

static val_ptr val_new_field(field_ptr f) {
  val_ptr v = pbc_malloc(sizeof(*v));
  v->type = v_field;
  v->data = f;
  return v;
}

static val_ptr fun_bin(
    void (*binop)(element_ptr, element_ptr, element_ptr),
    tree_ptr t) {
  // TODO: Check there are two args, types match.
  val_ptr v0 = tree_eval(darray_at(t->child, 0));
  val_ptr v1 = tree_eval(darray_at(t->child, 1));
  binop(v0->data, v0->data, v1->data);
  return v0;
}

val_ptr fun_add(tree_ptr t) { return fun_bin(element_add, t); }
val_ptr fun_sub(tree_ptr t) { return fun_bin(element_sub, t); }
val_ptr fun_mul(tree_ptr t) { return fun_bin(element_mul, t); }
val_ptr fun_div(tree_ptr t) { return fun_bin(element_div, t); }
val_ptr fun_pow(tree_ptr t) { return fun_bin(element_pow_zn, t); }

val_ptr fun_self(tree_ptr t) {
  // TODO: Write element_clone(), or at least element_new().
  element_ptr e = pbc_malloc(sizeof(*e));
  element_init_same_as(e, t->data);
  element_set(e, t->data);
  return val_new_element(e);
}

val_ptr fun_list(tree_ptr t) {
  element_ptr e = NULL;
  int n = darray_count(t->child);
  int i;
  for(i = 0; i < n; i++) {
    val_ptr x = tree_eval(darray_at(t->child, i));
    // TODO: Check x is element.
    if (!i) e = multiz_new_list(x->data);
    else multiz_append(e, x->data);
  }
  return val_new_element(e);
}

static val_ptr fun_id(tree_ptr t) {
  val_ptr x = symtab_at(reserved, t->data);
  if (!x) x = symtab_at(tab, t->data);
  if (!x) return NULL;  // TODO: Return error.
  return x->type->eval(x);
}

static val_ptr fun_fun(tree_ptr t) {
  val_ptr x = tree_eval(t->data);
  return x->type->funcall(x, t);
}

val_ptr fun_uminus(tree_ptr t) {
  val_ptr v = tree_eval(t->data);
  // TODO: Check v is an element.
  element_neg(v->data, v->data);
  return v;
}

val_ptr fun_assign(tree_ptr t) {
  val_ptr v = tree_eval(darray_at(t->child, 0));
  // TODO: Check ID is not reserved.
  symtab_put(tab, v, t->data);
  return v;
}

void assign_field(field_ptr f, const char* s) {
  symtab_put(tab, val_new_field(f), s);
}

void tree_append_multiz(tree_ptr t, tree_ptr m) {
  darray_append(t->child, m);
}

tree_ptr tree_new(val_ptr (*fun)(tree_ptr), void *data) {
  tree_ptr res = pbc_malloc(sizeof(*res));
  res->fun = fun;
  res->data = data;
  res->child = NULL;
  return res;
}

tree_ptr tree_new_z(const char* s) {
  element_ptr e = pbc_malloc(sizeof(*e));
  element_init(e, Z);
  element_set_str(e, s, 0);
  return tree_new(fun_self, e);
}

tree_ptr tree_new_list(tree_ptr first) {
  tree_ptr t = tree_new(fun_list, NULL);
  t->child = darray_new();
  darray_append(t->child, first);
  return t;
}

tree_ptr tree_new_id(const char* s) {
  return tree_new(fun_id, pbc_strdup(s));
}

tree_ptr tree_new_funcall(void) {
  tree_ptr t = tree_new(fun_fun, NULL);
  t->child = darray_new();
  return t;
}

tree_ptr tree_new_uminus(tree_ptr x) {
  tree_ptr t = tree_new(fun_uminus, x);
  return t;
}

void tree_set_fun(tree_ptr t, tree_ptr src) {
  free(t->data);
  t->data = src;
}

void tree_fun_append(tree_ptr f, tree_ptr p) {
  darray_append(f->child, p);
}

tree_ptr tree_new_bin(val_ptr (*fun)(tree_ptr), tree_ptr x, tree_ptr y) {
  tree_ptr t = tree_new(fun, NULL);
  t->child = darray_new();
  darray_append(t->child, x);
  darray_append(t->child, y);
  return t;
}

tree_ptr tree_new_assign(tree_ptr l, tree_ptr r) {
  // TODO: Check l's type.
  tree_ptr t = tree_new(fun_assign, l->data);
  t->child = darray_new();
  darray_append(t->child, r);
  return t;
}

// Top-level evaluation of a syntax tree node.
void tree_eval_stmt(tree_ptr t) {
  if (!t) return;
  val_ptr v = tree_eval(t);
  if (t->fun != fun_assign && v) {
    v->type->out_str(stdout, v);
    putchar('\n');
  }
}

static val_ptr fun_nextprime(tree_ptr t) {
  // TODO: Check args, x is an element.
  val_ptr x = tree_eval(darray_at(t->child, 0));
  element_ptr e = x->data;
  mpz_t z;
  mpz_init(z);
  element_to_mpz(z, e);
  mpz_nextprime(z, z);
  element_set_mpz(e, z);
  return x;
}

static val_ptr fun_order(tree_ptr t) {
  // TODO: Check args, x is a field.
  val_ptr x = tree_eval(darray_at(t->child, 0));
  val_ptr v = pbc_malloc(sizeof(*v));
  field_ptr f = x->data;
  element_ptr e = pbc_malloc(sizeof(*e));
  element_init(e, Z);
  element_set_mpz(e, f->order);
  v->type = v_elem;
  v->data = e;
  return v;
}

static val_ptr fun_random(tree_ptr t) {
  // TODO: Check args, x is a field.
  val_ptr x = tree_eval(darray_at(t->child, 0));
  val_ptr v = pbc_malloc(sizeof(*v));
  element_ptr e = pbc_malloc(sizeof(*e));
  element_init(e, x->data);
  element_random(e);
  v->type = v_elem;
  v->data = e;
  return v;
}

static val_ptr fun_sqrt(tree_ptr t) {
  // TODO: Check args, x is element, x is square.
  val_ptr x = tree_eval(darray_at(t->child, 0));
  element_sqrt(x->data, x->data);
  return x;
}

static val_ptr fun_inv(tree_ptr t) {
  // TODO: Check args, x is element, x is invertible.
  val_ptr x = tree_eval(darray_at(t->child, 0));
  element_invert(x->data, x->data);
  return x;
}

static val_ptr fun_type(tree_ptr t) {
  // TODO: Check args.
  val_ptr x = tree_eval(darray_at(t->child, 0));
  puts(x->type->name);
  return x;
}

static val_ptr fun_pairing(tree_ptr t) {
  // TODO: Check args.
  val_ptr x = tree_eval(darray_at(t->child, 0));
  val_ptr y = tree_eval(darray_at(t->child, 1));
  element_ptr xe = x->data;
  element_ptr e = element_new(xe->field->pairing->GT);
  bilinear_map(e, xe, y->data);
  return val_new_element(e);
}

static val_ptr fun_gf(tree_ptr t) {
  // TODO: Check args, x is an element.
  val_ptr x = tree_eval(darray_at(t->child, 0));
  element_ptr e = x->data;
  mpz_t z;
  mpz_init(z);
  element_to_mpz(z, e);
  field_ptr f = pbc_malloc(sizeof(*f));
  field_init_fp(f, z);
  mpz_clear(z);
  return val_new_field(f);
}

static void init_pairing(const char *s) {
  pairing_init_set_str(pairing, s);
  assign_field(pairing->G1, "G1");
  assign_field(pairing->G2, "G2");
  assign_field(pairing->GT, "GT");
  assign_field(pairing->Zr, "Zr");
}

static char *aparam =
"type a\n"
"q 8780710799663312522437781984754049815806883199414208211028653399266475630880222957078625179422662221423155858769582317459277713367317481324925129998224791\n"
"h 12016012264891146079388821366740534204802954401251311822919615131047207289359704531102844802183906537786776\n"
"r 730750818665451621361119245571504901405976559617\n"
"exp2 159\n"
"exp1 107\n"
"sign1 1\n"
"sign0 1\n";

static char *dparam =
"type d\n"
"q 625852803282871856053922297323874661378036491717\n"
"n 625852803282871856053923088432465995634661283063\n"
"h 3\n"
"r 208617601094290618684641029477488665211553761021\n"
"a 581595782028432961150765424293919699975513269268\n"
"b 517921465817243828776542439081147840953753552322\n"
"k 6\n"
"nk 60094290356408407130984161127310078516360031868417968262992864809623507269833854678414046779817844853757026858774966331434198257512457993293271849043664655146443229029069463392046837830267994222789160047337432075266619082657640364986415435746294498140589844832666082434658532589211525696\n"
"hk 1380801711862212484403205699005242141541629761433899149236405232528956996854655261075303661691995273080620762287276051361446528504633283152278831183711301329765591450680250000592437612973269056\n"
"coeff0 472731500571015189154958232321864199355792223347\n"
"coeff1 352243926696145937581894994871017455453604730246\n"
"coeff2 289113341693870057212775990719504267185772707305\n"
"nqr 431211441436589568382088865288592347194866189652\n";

static char *eparam =
"type e\n"
"q 7245986106510086080714203333362098431608853335867425877960916928496629182991629664903654100214900946450053872786629995869445693724001299041657434948257845644905153122838458864000479326695430719258600053239930483226650953770354174712511646273516974069245462534034085895319225452125649979474047163305307830001\n"
"r 730750862221594424981965739670091261094297337857\n"
"h 13569343110918781839835249021482970252603216587988030044836106948825516930173270978617489032334001006615524543925753725725046733884363846960470444404747241287743773746682188521738728797153760275116924829183670000\n"
"a 7130970454025799000067946137594446075551569949583815943390108723282396973737794273397246892274981883807989525599540630855644968426794929215599380425269625872763801485968007136000471718335185787206876242871042697778608875139078711621836858237429403052273312335081163896980825048123655535355411494046493419999\n"
"b 7169309004853894693616698536183663527570664411678352588247044791687141043489072737232715961588288238022010974661903752526911876859197052490952065266265699130144252031591491045333807587788600764557450846327338626261289568016170532652061787582791926724597362401398804563093625182790987016728290050466098223333\n"
"exp2 159\n"
"exp1 135\n"
"sign1 1\n"
"sign0 1\n";

static char *fparam =
"type f\n"
"q 205523667896953300194896352429254920972540065223\n"
"r 205523667896953300194895899082072403858390252929\n"
"b 40218105156867728698573668525883168222119515413\n"
"beta 115334401956802802075595682801335644058796914268\n"
"alpha0 191079354656274778837764015557338301375963168470\n"
"alpha1 71445317903696340296199556072836940741717506375\n";

static char *gparam =
"type g\n"
"q 503189899097385532598615948567975432740967203\n"
"n 503189899097385532598571084778608176410973351\n"
"h 1\n"
"r 503189899097385532598571084778608176410973351\n"
"a 465197998498440909244782433627180757481058321\n"
"b 463074517126110479409374670871346701448503064\n"
"k 10\n"
"nk 1040684643531490707494989587381629956832530311976146077888095795458709511789670022388326295177424065807612879371896982185473788988016190582073591316127396374860265835641044035656044524481121528846249501655527462202999638159773731830375673076317719519977183373353791119388388468745670818193868532404392452816602538968163226713846951514831917487400267590451867746120591750902040267826351982737642689423713163967384383105678367875981348397359466338807\n"
"hk 4110127713690841149713310614420858884651261781185442551927080083178682965171097172366598236129731931693425629387502221804555636704708008882811353539555915064049685663790355716130262332064327767695339422323460458479884756000782939428852120522712008037615051139080628734566850259704397643028017435446110322024094259858170303605703280329322675124728639532674407\n"
"coeff0 67343110967802947677845897216565803152319250\n"
"coeff1 115936772834120270862756636148166314916823221\n"
"coeff2 87387877425076080433559927080662339215696505\n"
"coeff3 433223145899090928132052677121692683015058909\n"
"coeff4 405367866213598664862417230702935310328613596\n"
"nqr 22204504160560785687198080413579021865783099\n";

static val_ptr fun_init_pairing_a(tree_ptr t) {
  UNUSED_VAR(t);
  init_pairing(aparam);
  return NULL;
}

static val_ptr fun_init_pairing_d(tree_ptr t) {
  UNUSED_VAR(t);
  init_pairing(dparam);
  return NULL;
}

static val_ptr fun_init_pairing_e(tree_ptr t) {
  UNUSED_VAR(t);
  init_pairing(eparam);
  return NULL;
}

static val_ptr fun_init_pairing_f(tree_ptr t) {
  UNUSED_VAR(t);
  init_pairing(fparam);
  return NULL;
}

static val_ptr fun_init_pairing_g(tree_ptr t) {
  UNUSED_VAR(t);
  init_pairing(gparam);
  return NULL;
}

static void builtin(val_ptr(*fun)(tree_ptr), const char *s) {
  val_ptr v = pbc_malloc(sizeof(*v));
  v->type = v_fun;
  v->data = fun;
  symtab_put(reserved, v, s);
}

int end_of_input;

int yywrap(void) {
  static char *currentline;
  static YY_BUFFER_STATE st;
  yy_delete_buffer(st);
  free(currentline);
  currentline = pbc_getline(option_prompt);
  if (!currentline) {
    end_of_input = 1;
    return 1;
  }
  int n = strlen(currentline);
  currentline = realloc(currentline, n + 2);
  currentline[n] = '\n';
  currentline[n + 1] = '\0';
  st = yy_scan_string(currentline);
  //if (option_echo) puts(currentline);
  return 0;
}

int main(int argc, char **argv) {
  for (;;) {
    int c = getopt(argc, argv, "y");
    if (c == -1) break;
    switch (c) {
      case 'y':
        option_easy = 1;
	option_prompt = "> ";
        break;
      default:
        fprintf(stderr, "unrecognized option: %c\n", c);
        break;
    }
  }

  field_init_multiz(Z);
  symtab_init(tab);

  builtin(fun_random, "rnd");
  builtin(fun_random, "random");
  builtin(fun_order, "ord");
  builtin(fun_order, "order");
  builtin(fun_nextprime, "nextprime");
  builtin(fun_sqrt, "sqrt");
  builtin(fun_inv, "inv");
  builtin(fun_type, "type");
  builtin(fun_pairing, "pairing");
  builtin(fun_gf, "GF");
  builtin(fun_init_pairing_a, "init_pairing_a");
  builtin(fun_init_pairing_d, "init_pairing_d");
  builtin(fun_init_pairing_e, "init_pairing_e");
  builtin(fun_init_pairing_f, "init_pairing_f");
  builtin(fun_init_pairing_g, "init_pairing_g");
  fun_init_pairing_a(NULL);

  symtab_put(reserved, val_new_field(Z), "Z");
  symtab_clear(tab);
  field_clear(Z);

  yywrap();
  while (!end_of_input) {
    if (2 == yyparse()) pbc_die("parser out of memory");
  }
  putchar('\n');
  return 0;
}
