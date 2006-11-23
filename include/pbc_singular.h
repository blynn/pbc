//requires
// * stdio.h
// * gmp.h
// * field.h
#ifndef __PBC_SINGULAR_H__
#define __PBC_SINGULAR_H__

void field_init_curve_singular_with_node(field_t c, field_t field);
void pairing_init_singular_with_node(pairing_t pairing, mpz_t q);

#endif //__PBC_SINGULAR_H__
