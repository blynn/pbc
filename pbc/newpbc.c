// Pairing-Based Calculator.

// TODO: Garbage collection.
#include "pbc.h"
#include "pbc_fieldmpz.h"
#include "pbc_fp.h"

#include "misc/darray.h"
#include "misc/symtab.h"

#include "pbc_tree.h"

#include "parser.tab.h"

static field_t Z;

struct tree_s {
  element_ptr (*fun)(tree_ptr);
  void *data;
  darray_t child;
};

void yyerror(char *s) {
  fprintf (stderr, "%s\n", s);
}

int yyparse(void);

static element_ptr fun_bin(
    void (*binop)(element_ptr, element_ptr, element_ptr),
    tree_ptr t) {
  // TODO: Check there are two args, types match.
  element_ptr e0 = tree_eval(darray_at(t->child, 0));
  element_ptr e1 = tree_eval(darray_at(t->child, 1));
  binop(e0, e0, e1);
  return e0;
}

element_ptr fun_add(tree_ptr t) { return fun_bin(element_add, t); }
element_ptr fun_sub(tree_ptr t) { return fun_bin(element_sub, t); }
element_ptr fun_mul(tree_ptr t) { return fun_bin(element_mul, t); }
element_ptr fun_div(tree_ptr t) { return fun_bin(element_div, t); }
element_ptr fun_pow(tree_ptr t) { return fun_bin(element_pow_zn, t); }

element_ptr fun_self(tree_ptr t) {
  element_ptr e = pbc_malloc(sizeof(*e));
  element_init_same_as(e, t->data);
  element_set(e, t->data);
  return e;
}

tree_ptr tree_new(element_ptr (*fun)(tree_ptr), void *data) {
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

tree_ptr tree_new_bin(element_ptr (*fun)(tree_ptr), tree_ptr x, tree_ptr y) {
  tree_ptr t = tree_new(fun, NULL);
  darray_append(t->child, x);
  darray_append(t->child, y);
  return t;
}

// Top-level evaluation of a syntax tree node.
void tree_eval_stmt(tree_ptr t) {
  element_ptr e = tree_eval(t);
  //if (t != fun_assign) {
  element_out_str(stdout, 0, e);
  putchar('\n');
}

element_ptr tree_eval(tree_ptr t) {
  return t->fun(t);
}

int main(int argc, char **argv) {
  field_init_z(Z);
  yyparse();
  field_clear(Z);
  return 0;
}
