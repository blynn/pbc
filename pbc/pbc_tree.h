// Requires:
// * field.h
struct tree_s;
typedef struct tree_s *tree_ptr;
tree_ptr tree_new_z(const char* s);
tree_ptr tree_new_list(tree_ptr t);
void tree_append_multiz(tree_ptr t, tree_ptr m);
tree_ptr tree_new_id(const char* s);
tree_ptr tree_new_assign(tree_ptr l, tree_ptr r);
tree_ptr tree_new_funcall(void);
void tree_fun_append(tree_ptr f, tree_ptr p);
void tree_set_fun(tree_ptr dst, tree_ptr src);
void tree_eval_stmt(tree_ptr t);

tree_ptr tree_new_neg(tree_ptr t);
tree_ptr tree_new_add(tree_ptr x, tree_ptr y);
tree_ptr tree_new_sub(tree_ptr x, tree_ptr y);
tree_ptr tree_new_mul(tree_ptr x, tree_ptr y);
tree_ptr tree_new_div(tree_ptr x, tree_ptr y);
tree_ptr tree_new_pow(tree_ptr x, tree_ptr y);
tree_ptr tree_new_eq(tree_ptr x, tree_ptr y);
tree_ptr tree_new_ne(tree_ptr x, tree_ptr y);
tree_ptr tree_new_le(tree_ptr x, tree_ptr y);
tree_ptr tree_new_ge(tree_ptr x, tree_ptr y);
tree_ptr tree_new_lt(tree_ptr x, tree_ptr y);
tree_ptr tree_new_gt(tree_ptr x, tree_ptr y);
tree_ptr tree_new_ternary(tree_ptr cond, tree_ptr t1, tree_ptr t2);
