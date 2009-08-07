// Pairing-Based Calculator.

// TODO: Garbage collection.
#include "pbc.h"
#include "pbc_fieldmpz.h"
#include "pbc_fp.h"

#include "misc/darray.h"
#include "misc/symtab.h"

#include "pbc_tree.h"

#include "parser.tab.h"

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

void yyerror(char *s) {
  fprintf (stderr, "%s\n", s);
}

int yyparse(void);

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
  element_ptr e = pbc_malloc(sizeof(*e));
  element_init_same_as(e, t->data);
  element_set(e, t->data);
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

tree_ptr tree_new(val_ptr (*fun)(tree_ptr), void *data) {
  tree_ptr res = pbc_malloc(sizeof(*res));
  res->fun = fun;
  res->data = data;
  res->child = NULL;
  return res;
}

tree_ptr tree_new_z(const char* s) {
  element_ptr z = pbc_malloc(sizeof(*z));
  element_init(z, Z);
  element_set_str(z, s, 0);
  return tree_new(fun_self, z);
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
  val_ptr v = tree_eval(t);
  if (t->fun != fun_assign) {
    v->type->out_str(stdout, v);
    putchar('\n');
  }
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

static void builtin(val_ptr(*fun)(tree_ptr), const char *s) {
  val_ptr v = pbc_malloc(sizeof(*v));
  v->type = v_fun;
  v->data = fun;
  symtab_put(reserved, v, s);
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

int main(int argc, char **argv) {
  field_init_z(Z);
  symtab_init(tab);

  builtin(fun_random, "rnd");
  builtin(fun_random, "random");
  builtin(fun_sqrt, "sqrt");
  builtin(fun_inv, "inv");
  builtin(fun_type, "type");
  builtin(fun_pairing, "pairing");
  pairing_init_set_str(pairing, aparam);
  assign_field(pairing->G1, "G1");
  assign_field(pairing->G2, "G2");
  assign_field(pairing->GT, "GT");
  assign_field(pairing->Zr, "Zr");

  symtab_put(reserved, val_new_field(Z), "Z");
  yyparse();
  symtab_clear(tab);
  field_clear(Z);
  return 0;
}
