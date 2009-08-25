/*
 * Quadratic field extensions.
 */

//requires
// * field.h
#ifndef __PBC_FIELDQUADRATIC_H__
#define __PBC_FIELDQUADRATIC_H__

// Initialize L as K[sqrt(a)], where a is a quadratic nonresidue of K. We
// automatically randomly generate a if necessary (see field_get_nqr() in
// field.c).
void field_init_quadratic(field_ptr L, field_ptr K);

// Initialize L as K[i], where i = sqrt(-1). Faster than the generic version.
// Requires -1 to be a quadratic nonresidue in K.
void field_init_fi(field_ptr L, field_ptr K);

// Naturally map an element from a field K to K[a].
void element_field_to_quadratic(element_ptr out, element_ptr in);
void element_field_to_fi(element_ptr a, element_ptr b);

#endif //__PBC_FIELDQUADRATIC_H__
