// some ternary extension fields,
// including $GF(3^m) = GF(3)[x]/(x^m + x^t + 2)$,
//           $GF(3^{2*m}) = GF(3^m)[x]/(x^2 + 1)$,
//           $GF(3^{3*m}) = GF(3^m)[x]/(x^3 - x - 1)$,
//           and $GF(3^{6*m}) = GF(3^{2*m})[x]/(x^3 - x - 1)$
//
// Requires:
// * pbc_field.h

#ifndef __PBC_TERNARY_EXTENSION_FIELD_H__
#define __PBC_TERNARY_EXTENSION_FIELD_H__

/* initialize $f$ as $GF(3)[x]/(x^m + x^t + 2)$ */
void field_init_gf3m(field_t f, unsigned m, unsigned t);

/* initialize $f$ as $base_field[x]/(x^2 + 1)$ */
void field_init_gf32m(field_t f, field_t base_field);

/* initialize $f$ as $base_field[x]/(x^3 - x - 1)$ */
void field_init_gf33m(field_t f, field_t base_field);

#endif //__PBC_TERNARY_EXTENSION_FIELD_H__
