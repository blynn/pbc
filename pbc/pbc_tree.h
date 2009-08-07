// Requires:
// * field.h
struct tree_s;
typedef struct tree_s *tree_ptr;
struct val_s;
typedef struct val_s *val_ptr;
tree_ptr tree_new_z(const char* s);
tree_ptr tree_new_id(const char* s);
tree_ptr tree_new_assign(tree_ptr l, tree_ptr r);
tree_ptr tree_new_uminus(tree_ptr x);
tree_ptr tree_new_bin(val_ptr (*)(tree_ptr), tree_ptr x, tree_ptr y);
tree_ptr tree_new_funcall(void);
void tree_fun_append(tree_ptr f, tree_ptr p);
void tree_set_fun(tree_ptr dst, tree_ptr src);
void tree_eval_stmt(tree_ptr t);
val_ptr fun_add(tree_ptr t);
val_ptr fun_sub(tree_ptr t);
val_ptr fun_mul(tree_ptr t);
val_ptr fun_div(tree_ptr t);
val_ptr fun_pow(tree_ptr t);
