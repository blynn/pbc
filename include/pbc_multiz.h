// Multinomnials with integer coefficients.

//requires
// * field.h

#ifndef __PBC_FIELDMULTI_H__
#define __PBC_FIELDMULTI_H__

void field_init_multiz(field_ptr f);

element_ptr multiz_new_list(element_ptr e);
void multiz_append(element_ptr l, element_ptr m);

#endif //__PBC_FIELDMULTI_H__
