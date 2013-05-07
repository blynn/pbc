// Pairing-Based Calculator.

// TODO: Garbage collection.
// TODO: Recursion (stack frames), anonymous functions.

#include <unistd.h>  // For getopt.

#include "pbc.h"
#include "pbc_fp.h"
#include "pbc_z.h"
#include "pbc_multiz.h"
#include "pbc_poly.h"

#include "misc/darray.h"
#include "misc/symtab.h"

#include "pbc_tree.h"

#include "lex.yy.h"
#include "parser.tab.h"

int option_easy = 0;
const char *option_prompt;

char *pbc_getline(const char *prompt);

void yyerror(char *s) { fprintf(stderr, "%s\n", s); }
int yyparse(void);

// Symbol table holding built-in functions and variables.
static symtab_t reserved;
// Symbol table holding user-defined variable and function names.
static symtab_t tab;

static field_t M;
static field_t Z;
static pairing_t pairing;

struct val_s;
typedef struct val_s *val_ptr;

struct fun_s;
typedef struct fun_s *fun_ptr;

// Syntax tree node.
struct tree_s {
  // Evaluates this node.
  val_ptr (*eval)(tree_ptr);
  union {
    const char *id;
    element_ptr elem;
    // Built-in function.
    fun_ptr fun;
    // Child nodes.
    darray_ptr child;
  };
};

enum {
  ARITY_VARIABLE = -1,
};

// The interface of a val_ptr shared amongst many val_ptr objects.
// Analog of C++ class.
struct val_type_s {
  // One of element, field, function, error.
  char *name;
  // Print out current value.
  void (*out_str)(FILE *, val_ptr);
  // Called when a variable is evaluated, e.g. "foo;".
  val_ptr (*eval)(val_ptr);
  // Called when a variable is used as a function, e.g. "foo();".
  val_ptr (*funcall)(val_ptr, tree_ptr);
};

// Functions plus type checking data.
struct fun_s {
  const char *name;
  val_ptr (*run)(val_ptr[]);
  int arity;
  const struct val_type_s **sig;
};
typedef struct fun_s fun_t[1];

// When interpreting, each node of the syntax tree recursively evaluates
// its children then returns a val_ptr.
struct val_s {
  struct val_type_s *type;
  union {
    element_ptr elem;
    // User-defined function.
    tree_ptr def;
    // Built-in function.
    fun_ptr fun;
    field_ptr field;
    const char *msg;
  };
};

static val_ptr val_new_element(element_ptr e);
static val_ptr val_new_field(field_ptr e);
static val_ptr val_new_error(const char *msg, ...);

// Evaluates syntax tree node.
static val_ptr tree_eval(tree_ptr t) {
  return t->eval(t);
}

static void v_elem_out(FILE* stream, val_ptr v) {
  element_out_str(stream, 0, v->elem);
}

static val_ptr v_elem_eval(val_ptr v) {
  element_ptr e = pbc_malloc(sizeof(*e));
  element_init_same_as(e, v->elem);
  element_set(e, v->elem);
  return val_new_element(e);
}

static void v_builtin_out(FILE* stream, val_ptr v) {
  // TODO: Print types of arguments.
  fprintf(stream, "built-in function %s, arity %d",
      v->fun->name, v->fun->arity);
}

static void v_define_out(FILE* stream, val_ptr v) {
  fprintf(stream, "user-defined function %s",
      ((tree_ptr) darray_at(v->def->child, 0))->id);
}

static val_ptr v_builtin(val_ptr v, tree_ptr t) {
  fun_ptr fun = v->fun;
  int n = fun->arity;
  if (1 + n != darray_count(t->child)) {
    return val_new_error("%s: wrong number of arguments", fun->name);
  }
  val_ptr arg[n];
  int i;
  for(i = 0; i < n; i++) {
    arg[i] = tree_eval(darray_at(t->child, i));
    if (fun->sig[i] && arg[i]->type != fun->sig[i]) {
      return val_new_error("%s: argument %d type mismatch", fun->name, i + 1);
    }
  }
  return fun->run(arg);
}

static void eval_stmt(void *ptr) {
  tree_eval(ptr);
}

static val_ptr v_def_call(val_ptr v, tree_ptr t) {
  int i;
  const char* name = ((tree_ptr) darray_at(v->def->child, 0))->id;
  darray_ptr parm = ((tree_ptr) darray_at(v->def->child, 1))->child;
  int n = darray_count(parm);
  if (1 + n != darray_count(t->child)) {
    return val_new_error("%s: wrong number of arguments", name);
  }
  for(i = 0; i < n; i++) {
    const char *id = ((tree_ptr) darray_at(parm, i))->id;
    val_ptr v1 = tree_eval(darray_at(t->child, i));
    // TODO: Stack frames for recursion.
    symtab_put(tab, v1, id);
  }
  // Evaluate function body.
  darray_ptr a = ((tree_ptr) darray_at(v->def->child, 2))->child;
  darray_forall(a, eval_stmt);
  return NULL;
}

static val_ptr v_field_cast(val_ptr v, tree_ptr t) {
  // TODO: Check args, x is an element.
  val_ptr x = tree_eval(darray_at(t->child, 0));
  element_ptr e = x->elem;
  if (e->field == M) {
    if (v->field == M) return x;
    element_ptr e2 = element_new(v->field);
    if (element_is0(e)) // if 'set0' is not 'set1' in base field of GT, but we hope 'GT(0)' calls 'set1', we may directly call 'element_set0' here
      element_set0(e2);
    else if (element_is1(e)) // reason is same as above
      element_set1(e2);
    else
      element_set_multiz(e2, e->data);
    x->elem = e2;
    return x;
  }
  if (v->field == M) {
    // Map to/from integer. TODO: Map to/from multiz instead.
    mpz_t z;
    mpz_init(z);
    element_to_mpz(z, e);
    element_clear(e);
    element_init(e, v->field);
    element_set_mpz(e, z);
    mpz_clear(z);
  }
  return x;
}

static void v_field_out(FILE* stream, val_ptr v) {
  field_out_info(stream, v->field);
}

static val_ptr v_self(val_ptr v) {
  return v;
}

static void v_err_out(FILE* stream, val_ptr v) {
  fprintf(stream, "%s", v->msg);
}

static val_ptr v_errcall(val_ptr v, tree_ptr t) {
  UNUSED_VAR(t);
  return v;
}

static struct val_type_s
  // TODO: Replace NULL with get_coeff.
  v_elem[1]  = {{  "element",    v_elem_out, v_elem_eval, NULL }},
  v_field[1] = {{    "field",   v_field_out,      v_self, v_field_cast }},
  v_fun[1]   = {{  "builtin", v_builtin_out,      v_self, v_builtin }},
  v_def[1]   = {{ "function",  v_define_out,      v_self, v_def_call }},
  v_error[1] = {{    "error",     v_err_out,      v_self, v_errcall }};

// Function signature constants for type checking.
const struct val_type_s *sig_field[] = { v_field };
const struct val_type_s *sig_elem[] = { v_elem };
const struct val_type_s *sig_any[] = { NULL };
const struct val_type_s *sig_elem_elem[] = { v_elem, v_elem };
const struct val_type_s *sig_field_elem[] = { v_field, v_elem };

static val_ptr val_new_element(element_ptr e) {
  val_ptr v = pbc_malloc(sizeof(*v));
  v->type = v_elem;
  v->elem = e;
  return v;
}

static val_ptr val_new_field(field_ptr f) {
  val_ptr v = pbc_malloc(sizeof(*v));
  v->type = v_field;
  v->field = f;
  return v;
}

static val_ptr val_new_error(const char *msg, ...) {
  va_list params;
  char buf[80];

  va_start(params, msg);
  vsnprintf(buf, 80, msg, params);
  va_end(params);

  val_ptr v = pbc_malloc(sizeof(*v));
  v->type = v_error;
  v->msg = pbc_strdup(buf);
  return v;
}

static val_ptr val_new_fun(fun_ptr fun) {
  val_ptr v = pbc_malloc(sizeof(*v));
  v->type = v_fun;
  v->fun = fun;
  return v;
}

static val_ptr fun_bin(
    void (*binop)(element_ptr, element_ptr, element_ptr),
    val_ptr v[]) {
  binop(v[0]->elem, v[0]->elem, v[1]->elem);
  return v[0];
}

static val_ptr run_add(val_ptr v[]) { return fun_bin(element_add, v); }
static val_ptr run_sub(val_ptr v[]) { return fun_bin(element_sub, v); }
static val_ptr run_mul(val_ptr v[]) { return fun_bin(element_mul, v); }
static val_ptr run_div(val_ptr v[]) { return fun_bin(element_div, v); }
static val_ptr run_pow(val_ptr v[]) { return fun_bin(element_pow_zn, v); }

static fun_t fun_add = {{ "add", run_add, 2, sig_elem_elem }};
static fun_t fun_sub = {{ "sub", run_sub, 2, sig_elem_elem }};
static fun_t fun_mul = {{ "mul", run_mul, 2, sig_elem_elem }};
static fun_t fun_div = {{ "div", run_div, 2, sig_elem_elem }};
static fun_t fun_pow = {{ "pow", run_pow, 2, sig_elem_elem }};

static val_ptr fun_cmp(val_ptr v[], int (*fun)(int)) {
  int i = element_cmp(v[0]->elem, v[1]->elem);
  element_ptr e = pbc_malloc(sizeof(*e));
  element_init(e, M);
  element_set_si(e, fun(i));
  v[0]->elem = e;
  return v[0];
}

static int is0(int i) {
  return i == 0;
}

static int isnot0(int i) {
  return i != 0;
}

static int isle(int i) {
  return i <= 0;
}

static int isge(int i) {
  return i >= 0;
}

static int islt(int i) {
  return i < 0;
}

static int isgt(int i) {
  return i > 0;
}

static val_ptr run_eq(val_ptr v[]) {
  return fun_cmp(v, is0);
}

static val_ptr run_ne(val_ptr v[]) {
  return fun_cmp(v, isnot0);
}

static val_ptr run_le(val_ptr v[]) {
  return fun_cmp(v, isle);
}

static val_ptr run_ge(val_ptr v[]) {
  return fun_cmp(v, isge);
}
static val_ptr run_lt(val_ptr v[]) {
  return fun_cmp(v, islt);
}
static val_ptr run_gt(val_ptr v[]) {
  return fun_cmp(v, isgt);
}

static fun_t fun_eq = {{ "==", run_eq, 2, sig_elem_elem }};
static fun_t fun_ne = {{ "!=", run_ne, 2, sig_elem_elem }};
static fun_t fun_le = {{ "<=", run_le, 2, sig_elem_elem }};
static fun_t fun_ge = {{ ">=", run_ge, 2, sig_elem_elem }};
static fun_t fun_lt = {{ "<", run_lt, 2, sig_elem_elem }};
static fun_t fun_gt = {{ ">", run_gt, 2, sig_elem_elem }};

static val_ptr eval_elem(tree_ptr t) {
  // TODO: Write element_clone(), or at least element_new().
  element_ptr e = pbc_malloc(sizeof(*e));
  element_init_same_as(e, t->elem);
  element_set(e, t->elem);
  return val_new_element(e);
}

static val_ptr eval_list(tree_ptr t) {
  element_ptr e = NULL;
  int n = darray_count(t->child);
  int i;
  for(i = 0; i < n; i++) {
    val_ptr x = tree_eval(darray_at(t->child, i));
    // TODO: Also check x is a multiz.
    if (v_error == x->type) {
      return x;
    }
    if (v_elem != x->type) {
      return val_new_error("element expected in list");
    }
    if (!i) e = multiz_new_list(x->elem);
    else multiz_append(e, x->elem);
  }
  return val_new_element(e);
}

static val_ptr eval_ternary(tree_ptr t) {
  val_ptr x = tree_eval(darray_at(t->child, 0));
  if (v_error == x->type) {
    return x;
  }
  if (x->type != v_elem) {
    return val_new_error("element expected in ternary operator");
  }
  if (!element_is0(x->elem)) {
    return tree_eval(darray_at(t->child, 1));
  }
  return tree_eval(darray_at(t->child, 2));
}

static val_ptr eval_id(tree_ptr t) {
  val_ptr x = symtab_at(reserved, t->id);
  if (!x) x = symtab_at(tab, t->id);
  if (!x) {
    return val_new_error("undefined variable %s", t->id);
  }
  return x->type->eval(x);
}

static val_ptr eval_funcall(tree_ptr t) {
  val_ptr x = tree_eval(darray_last(t->child));
  return x->type->funcall(x, t);
}

static val_ptr eval_fun(tree_ptr t) {
  return val_new_fun(t->fun);
}

static val_ptr run_neg(val_ptr v[]) {
  element_neg(v[0]->elem, v[0]->elem);
  return v[0];
}
static fun_t fun_neg = {{ "neg", run_neg, 1, sig_elem }};

static val_ptr eval_assign(tree_ptr t) {
  tree_ptr tid = darray_at(t->child, 0);
  val_ptr v = tree_eval(darray_at(t->child, 1));
  if (symtab_at(reserved, tid->id)) {
    return val_new_error("%s is reserved", tid->id);
  }
  symtab_put(tab, v, tid->id);
  return v;
}

static void assign_field(field_ptr f, const char* s) {
  symtab_put(tab, val_new_field(f), s);
}

tree_ptr tree_new(val_ptr (*eval)(tree_ptr)) {
  tree_ptr res = pbc_malloc(sizeof(*res));
  res->eval = eval;
  return res;
}

tree_ptr tree_new_z(const char* s) {
  element_ptr e = pbc_malloc(sizeof(*e));
  element_init(e, M);
  element_set_str(e, s, 0);
  tree_ptr t = tree_new(eval_elem);
  t->elem = e;
  return t;
}

static val_ptr eval_err(tree_ptr t) {
  UNUSED_VAR(t);
  pbc_die("BUG: shouldn't reach here!");
}

tree_ptr tree_new_empty_stmt_list() {
  tree_ptr t = tree_new(eval_err);
  t->child = darray_new();
  return t;
}

tree_ptr tree_new_empty_parms() {
  tree_ptr t = tree_new(eval_err);
  t->child = darray_new();
  return t;
}

static val_ptr eval_define(tree_ptr t) {
  val_ptr v = pbc_malloc(sizeof(*v));
  v->type = v_def;
  v->def = t;
  symtab_put(tab, v, ((tree_ptr) darray_at(t->child, 0))->id);
  return v;
}

tree_ptr tree_new_define(tree_ptr id, tree_ptr parm, tree_ptr body) {
  tree_ptr t = tree_new(eval_define);
  t->child = darray_new();
  darray_append(t->child, id);
  darray_append(t->child, parm);
  darray_append(t->child, body);
  return t;
}

tree_ptr tree_new_list(tree_ptr first) {
  tree_ptr t = tree_new(eval_list);
  t->child = darray_new();
  darray_append(t->child, first);
  return t;
}

tree_ptr tree_new_ternary(tree_ptr cond, tree_ptr t1, tree_ptr t2) {
  tree_ptr t = tree_new(eval_ternary);
  t->child = darray_new();
  darray_append(t->child, cond);
  darray_append(t->child, t1);
  darray_append(t->child, t2);
  return t;
}

tree_ptr tree_new_id(const char* s) {
  tree_ptr t = tree_new(eval_id);
  t->id = pbc_strdup(s);
  return t;
}

tree_ptr tree_new_funcall(void) {
  tree_ptr t = tree_new(eval_funcall);
  t->child = darray_new();
  return t;
}

static tree_ptr tree_new_fun(fun_ptr fun) {
  tree_ptr t = tree_new(eval_fun);
  t->fun = fun;
  return t;
}

void tree_set_fun(tree_ptr f, tree_ptr src) {
  darray_append(f->child, src);
}

void tree_append(tree_ptr f, tree_ptr p) {
  darray_append(f->child, p);
}

tree_ptr tree_new_binary(fun_ptr fun, tree_ptr x, tree_ptr y) {
  tree_ptr t = tree_new_funcall();
  tree_append(t, x);
  tree_append(t, y);
  tree_set_fun(t, tree_new_fun(fun));
  return t;
}

static tree_ptr tree_new_unary(fun_ptr fun, tree_ptr x) {
  tree_ptr t = tree_new_funcall();
  tree_append(t, x);
  tree_set_fun(t, tree_new_fun(fun));
  return t;
}

tree_ptr tree_new_neg(tree_ptr t) {
  return tree_new_unary(fun_neg, t);
}
tree_ptr tree_new_add(tree_ptr x, tree_ptr y) {
  return tree_new_binary(fun_add, x, y);
}
tree_ptr tree_new_sub(tree_ptr x, tree_ptr y) {
  return tree_new_binary(fun_sub, x, y);
}
tree_ptr tree_new_mul(tree_ptr x, tree_ptr y) {
  return tree_new_binary(fun_mul, x, y);
}
tree_ptr tree_new_div(tree_ptr x, tree_ptr y) {
  return tree_new_binary(fun_div, x, y);
}
tree_ptr tree_new_pow(tree_ptr x, tree_ptr y) {
  return tree_new_binary(fun_pow, x, y);
}
tree_ptr tree_new_eq(tree_ptr x, tree_ptr y) {
  return tree_new_binary(fun_eq, x, y);
}
tree_ptr tree_new_ne(tree_ptr x, tree_ptr y) {
  return tree_new_binary(fun_ne, x, y);
}
tree_ptr tree_new_le(tree_ptr x, tree_ptr y) {
  return tree_new_binary(fun_le, x, y);
}
tree_ptr tree_new_ge(tree_ptr x, tree_ptr y) {
  return tree_new_binary(fun_ge, x, y);
}
tree_ptr tree_new_lt(tree_ptr x, tree_ptr y) {
  return tree_new_binary(fun_lt, x, y);
}
tree_ptr tree_new_gt(tree_ptr x, tree_ptr y) {
  return tree_new_binary(fun_gt, x, y);
}

static val_ptr run_item(val_ptr v[]) {
  mpz_t z;
  mpz_init(z);
  element_to_mpz(z, v[1]->elem);
  int i = mpz_get_si(z);
  mpz_clear(z);
  element_ptr a = element_item(v[0]->elem, i);
  element_ptr e = pbc_malloc(sizeof(*e));
  element_init_same_as(e, a);
  element_set(e, a);
  return val_new_element(e);
}
static fun_t fun_item = {{ "item", run_item, 2, sig_elem_elem }};
tree_ptr tree_new_item(tree_ptr x, tree_ptr y) {
  return tree_new_binary(fun_item, x, y);
}

tree_ptr tree_new_assign(tree_ptr l, tree_ptr r) {
  // TODO: Check l's type.
  tree_ptr t = tree_new(eval_assign);
  t->child = darray_new();
  darray_append(t->child, l);
  darray_append(t->child, r);
  return t;
}

// Evaluate statement.
void tree_eval_stmt(tree_ptr stmt) {
  val_ptr v = tree_eval(stmt);
  if (v && v_error == v->type) {
    v->type->out_str(stdout, v);
    putchar('\n');
  } else if (stmt->eval != eval_assign && v) {
    v->type->out_str(stdout, v);
    putchar('\n');
  }
}

static val_ptr run_nextprime(val_ptr v[]) {
  element_ptr e = v[0]->elem;
  mpz_t z;
  mpz_init(z);
  element_to_mpz(z, e);
  mpz_nextprime(z, z);
  element_set_mpz(e, z);
  return v[0];
}
static fun_t fun_nextprime = {{ "nextprime", run_nextprime, 1, sig_elem }};

static val_ptr run_order(val_ptr v[]) {
  field_ptr f = v[0]->field;
  element_ptr e = pbc_malloc(sizeof(*e));
  element_init(e, M);
  element_set_mpz(e, f->order);
  return val_new_element(e);
}
static fun_t fun_ord = {{ "ord", run_order, 1, sig_field }};
static fun_t fun_order = {{ "order", run_order, 1, sig_field }};

static val_ptr run_random(val_ptr v[]) {
  element_ptr e = pbc_malloc(sizeof(*e));
  element_init(e, v[0]->field);
  element_random(e);
  return val_new_element(e);
}
static fun_t fun_rnd = {{ "rnd", run_random, 1, sig_field }};
static fun_t fun_random = {{ "random", run_random, 1, sig_field }};

static val_ptr run_sqrt(val_ptr v[]) {
  // TODO: Check v[0] is square.
  element_sqrt(v[0]->elem, v[0]->elem);
  return v[0];
}
static fun_t fun_sqrt = {{ "sqrt", run_sqrt, 1, sig_elem }};

static val_ptr run_invert(val_ptr v[]) {
  // TODO: Check v[0] is invertible.
  element_invert(v[0]->elem, v[0]->elem);
  return v[0];
}
static fun_t fun_inv = {{ "inv", run_invert, 1, sig_elem }};

static val_ptr run_type(val_ptr v[]) {
  puts(v[0]->type->name);
  return v[0];
}
static fun_t fun_type = {{ "type", run_type, 1, sig_any }};

static val_ptr run_pairing(val_ptr v[]) {
  element_ptr x = v[0]->elem;
  element_ptr e = element_new(x->field->pairing->GT);
  element_pairing(e, x, v[1]->elem);
  return val_new_element(e);
}
static fun_t fun_pairing = {{ "pairing", run_pairing, 2, sig_elem_elem }};

static val_ptr run_zmod(val_ptr v[]) {
  element_ptr e = v[0]->elem;
  mpz_t z;
  mpz_init(z);
  element_to_mpz(z, e);
  field_ptr f = pbc_malloc(sizeof(*f));
  field_init_fp(f, z);
  mpz_clear(z);
  return val_new_field(f);
}
static fun_t fun_zmod = {{ "zmod", run_zmod, 1, sig_elem }};

static val_ptr run_poly(val_ptr v[]) {
  field_ptr f = pbc_malloc(sizeof(*f));
  field_init_poly(f, v[0]->field);
  return val_new_field(f);
}
static fun_t fun_poly = {{ "poly", run_poly, 1, sig_field }};

static val_ptr run_polymod(val_ptr v[]) {
  // TODO: Check v[0] is a poly.
  field_ptr f = pbc_malloc(sizeof(*f));
  field_init_polymod(f, v[0]->elem);
  return val_new_field(f);
}
static fun_t fun_polymod = {{ "polymod", run_polymod, 1, sig_elem }};

static val_ptr run_extend(val_ptr v[]) {
  // TODO: Check v[1] is multiz poly.
  field_ptr fx = pbc_malloc(sizeof(*fx));
  field_init_poly(fx, v[0]->field);
  element_ptr poly = element_new(fx);
  element_set_multiz(poly, v[1]->elem->data);
  field_ptr f = pbc_malloc(sizeof(*f));
  field_init_polymod(f, poly);
  element_free(poly);
  return val_new_field(f);
}
static fun_t fun_extend = {{ "extend", run_extend, 1, sig_field_elem }};

static void init_pairing(const char *s) {
  pairing_init_set_str(pairing, s);
  assign_field(pairing->G1, "G1");
  assign_field(pairing->G2, "G2");
  assign_field(pairing->GT, "GT");
  assign_field(pairing->Zr, "Zr");
}

static val_ptr run_exit(val_ptr v[]) {
  mpz_t z;
  mpz_init(z);
  element_to_mpz(z, v[0]->elem);
  exit(mpz_get_si(z));
}
static fun_t fun_exit = {{ "exit", run_exit, 1, sig_elem }};

static val_ptr run_CHECK(val_ptr v[]) {
  if (element_is0(v[0]->elem)) {
    pbc_die("CHECK failed");
  }
  return v[0];
}
static fun_t fun_CHECK = {{ "CHECK", run_CHECK, 1, sig_elem }};

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

static char *iparam =
"type i\n"
"m 97\n"
"t 12\n"
"n 2726865189058261010774960798134976187171462721\n"
"n2 7\n";

static val_ptr run_init_pairing_a(val_ptr v[]) {
  UNUSED_VAR(v);
  init_pairing(aparam);
  return NULL;
}
static fun_t fun_init_pairing_a = {{
    "init_pairing_a", run_init_pairing_a, 0, NULL
    }};

static val_ptr run_init_pairing_d(val_ptr v[]) {
  UNUSED_VAR(v);
  init_pairing(dparam);
  return NULL;
}
static fun_t fun_init_pairing_d = {{
    "init_pairing_d", run_init_pairing_d, 0, NULL
    }};

static val_ptr run_init_pairing_e(val_ptr v[]) {
  UNUSED_VAR(v);
  init_pairing(eparam);
  return NULL;
}
static fun_t fun_init_pairing_e = {{
    "init_pairing_e", run_init_pairing_e, 0, NULL
    }};

static val_ptr run_init_pairing_f(val_ptr v[]) {
  UNUSED_VAR(v);
  init_pairing(fparam);
  return NULL;
}
static fun_t fun_init_pairing_f = {{
    "init_pairing_f", run_init_pairing_f, 0, NULL
    }};

static val_ptr run_init_pairing_g(val_ptr v[]) {
  UNUSED_VAR(v);
  init_pairing(gparam);
  return NULL;
}
static fun_t fun_init_pairing_g = {{
    "init_pairing_g", run_init_pairing_g, 0, NULL
    }};

static val_ptr run_init_pairing_i(val_ptr v[]) {
  UNUSED_VAR(v);
  init_pairing(iparam);
  return NULL;
}
static fun_t fun_init_pairing_i = {{
    "init_pairing_i", run_init_pairing_i, 0, NULL
    }};

static void builtin(fun_ptr fun) {
  symtab_put(reserved, val_new_fun(fun), fun->name);
}

int end_of_input;

int yywrap_return1(void) { return 1; }

int yywrap_readline(void) {
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

static int (*yywrapfun)(void);
int yywrap(void) {
  return yywrapfun();
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

  field_init_z(Z);
  field_init_multiz(M);
  symtab_init(tab);

  builtin(fun_rnd);
  builtin(fun_random);
  builtin(fun_ord);
  builtin(fun_order);
  builtin(fun_nextprime);
  builtin(fun_sqrt);
  builtin(fun_inv);
  builtin(fun_type);
  builtin(fun_pairing);
  builtin(fun_zmod);
  builtin(fun_poly);
  builtin(fun_polymod);
  builtin(fun_extend);
  builtin(fun_exit);
  builtin(fun_CHECK);
  builtin(fun_init_pairing_a);
  builtin(fun_init_pairing_d);
  builtin(fun_init_pairing_e);
  builtin(fun_init_pairing_f);
  builtin(fun_init_pairing_g);
  builtin(fun_init_pairing_i);
  run_init_pairing_a(NULL);
  symtab_put(reserved, val_new_field(M), "M");
  symtab_put(reserved, val_new_field(Z), "Z");

  if (argc > optind) {
    FILE *fp = fopen(argv[optind], "r");
    if (!fp) pbc_die("fopen failed on %s", argv[optind]);
    YY_BUFFER_STATE st = yy_create_buffer(fp, YY_BUF_SIZE);
    yy_switch_to_buffer(st);
    yywrapfun = yywrap_return1;
    yyparse();
    yy_delete_buffer(st);
  } else {
    yywrapfun = yywrap_readline;
    yywrap();
    while (!end_of_input) {
      if (2 == yyparse()) pbc_die("parser out of memory");
    }
    putchar('\n');
  }

  symtab_clear(tab);
  field_clear(M);
  return 0;
}
