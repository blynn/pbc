// Multinomnials with integer coefficients.

//requires
// * field.h

#ifndef __PBC_FIELDMULTI_H__
#define __PBC_FIELDMULTI_H__

void field_init_multiz(field_ptr f);

element_ptr multiz_new_list(element_ptr e);
void multiz_append(element_ptr l, element_ptr m);

void multiz_to_mpz(mpz_ptr z, multiz ep);
int multiz_is_z(multiz m);
multiz multiz_at(multiz m, int i);
int multiz_count(multiz m);
int multiz_is0(multiz m);

#endif //__PBC_FIELDMULTI_H__
