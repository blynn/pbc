// A subgroup of the group of points on an elliptic curve.
// Also used to represent quotient groups.
//
// We use the field_t structure even though E(K) is a group. Addition and
// multiplication both refer to the group operation.

// Requires:
// * stdio.h
// * gmp.h
// * field.h
#ifndef __PBC_CURVE_H__
#define __PBC_CURVE_H__

// Some initialization functions take an order parameter. This is meant to
// be the order of the subgroup, but might actually be the order of the twist.
// Certain routines initialize a curve, test a random point to see if it has
// the correct order, and if not, immediately twist the curve so that it does.
// TODO: Move such code into curve.c, so 'order' is always accurate.

// If cofac != NULL, then the field_t represents the subgroup of
// order = #E(K) / cofac.
//
// If not, and order = #E(K) then the field_t represents the entire E(K).
//
// Otherwise, if order is a factor of #E(K), then the field_t represents
// the quotient group of that order, namely E(K)/(#E(K)/order). No attempt is
// made to standardize the coset representative. This mode is useful for the
// Tate pairing (see thesis), where any coset representative of G2 suffices
// during the pairing computation.

// Initialize a subgroup of points on the curve Y^2 = X^3 + b.
void field_init_curve_b(field_ptr f, element_ptr b, mpz_t order, mpz_t cofac);

// Initialize a subgroup of points on the curve with the given j-invariant.
void field_init_curve_j(field_t f, element_ptr j, mpz_t order, mpz_t cofac);

// Initialize a subgroup of points on the curve Y^2 = X^3 + a X + b.
void field_init_curve_ab(field_ptr f, element_ptr a, element_ptr b, mpz_t order, mpz_t cofac);

// Reinitialize as the subgroup of points on the twist curve.
// Requires j-invariant of the original curve != 0, 1728.
// Mangles f, thus existing points of f become invalid.
// TODO: Refactor so we can remove this from the interface.
void field_reinit_curve_twist(field_t f);

// Compute trace of Frobenius at q^n given trace at q.
void pbc_mpz_trace_n(mpz_t res, mpz_t q, mpz_t trace, int n);

// Given q, t such that #E(F_q) = q - t + 1, compute #E(F_q^k).
void pbc_mpz_curve_order_extn(mpz_t res, mpz_t q, mpz_t t, int k);

void field_init_curve_with_map(field_ptr cnew, field_ptr c,
  field_ptr dstfield, fieldmap map);

void field_init_curve_ab_map(field_t cnew, field_t c,
  fieldmap map, field_ptr mapdest,
  mpz_t ordernew, mpz_t cofacnew);

void field_curve_use_random_solvefory(field_ptr f);

void field_curve_set_quotient_cmp(field_ptr c, mpz_t quotient_cmp);

#pragma GCC visibility push(hidden)
// Internal:

element_ptr curve_x_coord(element_t e);
element_ptr curve_y_coord(element_t e);
element_ptr curve_a_coeff(element_t e);
element_ptr curve_b_coeff(element_t e);
element_ptr curve_field_a_coeff(field_t f);
element_ptr curve_field_b_coeff(field_t f);

void curve_from_x(element_ptr e, element_t x);
void curve_set_si(element_t R, long int x, long int y);
void curve_set_gen_no_cofac(element_ptr a);

#pragma GCC visibility pop

#endif //__PBC_CURVE_H__
