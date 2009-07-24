// Polynomial rings R[x], and polynomial rings modulo polynomials,
// i.e. R[x]_{f(x)}.

// Requires:
// * gmp.h
// * field.h
// * darray.h
#ifndef __PBC_POLY_H__
#define __PBC_POLY_H__

int poly_coeff_count(element_ptr e);
element_ptr poly_coeff(element_ptr e, int i);

static inline int poly_degree(element_ptr e) {
  return poly_coeff_count(e) - 1;
}

field_ptr poly_base_field(element_t f);

void poly_alloc(element_ptr e, int n);
void poly_remove_leading_zeroes(element_ptr e);
void poly_set_coeff(element_ptr e, element_ptr a, int n);
void poly_setx(element_ptr f);
void poly_const_mul(element_ptr res, element_ptr a, element_ptr poly);
void poly_div(element_ptr quot, element_ptr rem,
    element_ptr a, element_ptr b);

void field_init_poly(field_ptr f, field_ptr base_field);
void field_init_polymod(field_ptr f, element_ptr poly);

void poly_gcd(element_ptr d, element_ptr f, element_ptr g);
// Return 1 if polynomial is irreducible, 0 otherwise.
int poly_is_irred(element_ptr f);
void poly_invert(element_ptr res, element_ptr f, element_ptr m);
void poly_random_monic(element_ptr f, int deg);
void poly_make_monic(element_t f, element_t g);
void element_field_to_poly(element_t poly, element_t constant);

void element_polymod_to_poly(element_ptr f, element_ptr e);
void element_field_to_polymod(element_ptr f, element_ptr a);
void element_poly_to_polymod_truncate(element_ptr f, element_ptr e);
element_ptr polymod_coeff(element_ptr e, int i);

void polymod_const_mul(element_ptr res, element_ptr a, element_ptr e);
int polymod_field_degree(field_t f);
#endif //__PBC_POLY_H__
