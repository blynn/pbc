// Polynomial rings R[x], and polynomial rings modulo polynomials,
// i.e. R[x]_{f(x)}.

// Requires:
// * gmp.h
// * field.h
// * darray.h
#ifndef __PBC_POLY_H__
#define __PBC_POLY_H__

// Returns deg f + 1.
int poly_coeff_count(element_ptr f);

// Returns coefficient of x^n in f.
// Assumes deg f >= n.
element_ptr poly_coeff(element_ptr f, int n);

// Returns deg f.
static inline int poly_degree(element_ptr f) {
  return poly_coeff_count(f) - 1;
}

// Returns base field of f (where the coefficients live).
field_ptr poly_base_field(element_t f);

void poly_alloc(element_ptr e, int n);

// Sets the coefficient of x^n to a.
void poly_set_coeff(element_ptr f, element_ptr a, int n);

// Sets f = x.
void poly_setx(element_ptr f);
void poly_const_mul(element_ptr res, element_ptr a, element_ptr poly);

// Initializes a polynomial ring.
void field_init_poly(field_ptr f, field_ptr base_field);

// Initializes a polynomial modulo ring.
// Requires poly to be monic.
void field_init_polymod(field_ptr f, element_ptr poly);

// TODO: Move findroot to poly.c and expose that instead of these functions:
// Sets d = gcd(f, g).
// Exposed because ecc/hilbert.c uses it in findroot.
void poly_gcd(element_ptr d, element_ptr f, element_ptr g);

// Sets f = c g for some constant c and is monic.
// Requires the leading coefficient of g to be a unit.
// Exposed because ecc/hilbert.c uses it in findroot.
void poly_make_monic(element_t f, element_t g);

// Returns 1 if polynomial is irreducible, 0 otherwise.
// Requires the polynomial to be monic.
int poly_is_irred(element_ptr f);
void poly_random_monic(element_ptr f, int deg);

void element_field_to_poly(element_t poly, element_t constant);

void element_polymod_to_poly(element_ptr f, element_ptr e);
void element_field_to_polymod(element_ptr f, element_ptr a);
void element_poly_to_polymod_truncate(element_ptr f, element_ptr e);
element_ptr polymod_coeff(element_ptr e, int i);

void polymod_const_mul(element_ptr res, element_ptr a, element_ptr e);
int polymod_field_degree(field_t f);
#endif //__PBC_POLY_H__
