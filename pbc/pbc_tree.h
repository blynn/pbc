// Requires:
// * field.h
struct tree_s;
typedef struct tree_s *tree_ptr;
tree_ptr tree_new_z(const char* s);
tree_ptr tree_new_bin(element_ptr (*)(tree_ptr), tree_ptr x, tree_ptr y);
void tree_eval_stmt(tree_ptr t);
element_ptr tree_eval(tree_ptr t);
element_ptr fun_add(tree_ptr t);
element_ptr fun_sub(tree_ptr t);
element_ptr fun_mul(tree_ptr t);
element_ptr fun_div(tree_ptr t);
element_ptr fun_pow(tree_ptr t);
