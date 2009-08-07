// Pairing-Based Calculator.

// TODO: Garbage collection.
#include "pbc.h"
#include "pbc_fieldmpz.h"
#include "pbc_fp.h"

#include "misc/darray.h"
#include "misc/symtab.h"

#include "pbc_tree.h"

#include "parser.tab.h"

// Main symbol table, holding variable and function names.
static symtab_t tab;

static field_t Z;

struct tree_s {
  val_ptr (*fun)(tree_ptr);
  void *data;
  darray_t child;
};

enum {
  V_ELEM,
  V_FIELD,
  V_ERROR,
};

struct val_s {
  int type;
  void *data;
};

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
  val_ptr v = pbc_malloc(sizeof(*v));
  element_ptr e = pbc_malloc(sizeof(*e));
  element_init_same_as(e, t->data);
  element_set(e, t->data);
  v->type = V_ELEM;
  v->data = e;
  return v;
}

val_ptr fun_id(tree_ptr t) {
  val_ptr x = symtab_at(tab, t->data);
  if (!x) return NULL;  // TODO: Return error.
  val_ptr v = pbc_malloc(sizeof(*v));
  element_ptr e = pbc_malloc(sizeof(*e));
  // TODO: Check x is V_ELEM.
  element_init_same_as(e, x->data);
  element_set(e, x->data);
  v->type = V_ELEM;
  v->data = e;
  return v;
}

val_ptr fun_assign(tree_ptr t) {
  val_ptr v = tree_eval(darray_at(t->child, 0));
  symtab_put(tab, v, t->data);
  return v;
}

tree_ptr tree_new(val_ptr (*fun)(tree_ptr), void *data) {
  tree_ptr res = pbc_malloc(sizeof(*res));
  res->fun = fun;
  res->data = data;
  darray_init(res->child);
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

tree_ptr tree_new_bin(val_ptr (*fun)(tree_ptr), tree_ptr x, tree_ptr y) {
  tree_ptr t = tree_new(fun, NULL);
  darray_append(t->child, x);
  darray_append(t->child, y);
  return t;
}

tree_ptr tree_new_assign(tree_ptr l, tree_ptr r) {
  // TODO: Check l's type.
  tree_ptr t = tree_new(fun_assign, l->data);
  darray_append(t->child, r);
  return t;
}

static void val_out_str(FILE* stream, val_ptr v) {
  switch(v->type) {
    case V_ELEM:
      element_out_str(stream, 0, v->data);
      break;
    default:
      pbc_warn("Unhandled V_ type");
      break;
  }
}

// Top-level evaluation of a syntax tree node.
void tree_eval_stmt(tree_ptr t) {
  val_ptr v = tree_eval(t);
  if (t->fun != fun_assign) {
    val_out_str(stdout, v);
    putchar('\n');
  }
}

val_ptr tree_eval(tree_ptr t) {
  return t->fun(t);
}

int main(int argc, char **argv) {
  field_init_z(Z);
  symtab_init(tab);
  yyparse();
  symtab_clear(tab);
  field_clear(Z);
  return 0;
}
