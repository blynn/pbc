//requires
// * stdio.h
// * gmp.h
// * field.h
#ifndef __PBC_CURVE_H__
#define __PBC_CURVE_H__

/* called in {e,f}_param.c */
void field_init_curve_b(field_ptr f, element_ptr b, mpz_t order, mpz_t cofac);

/* these are called in mnt.c */
void field_init_curve_j(field_t f, element_ptr j, mpz_t order, mpz_t cofac);

    //assumes j != 0, 1728
void twist_curve(field_ptr c);
    //WARNING: existing points will no longer lie on c
    //as this modifies c in place
void compute_trace_n(mpz_t res, mpz_t q, mpz_t trace, int n);
    //compute trace of Frobenius at q^n given trace at q
    //see p.105 of Blake, Seroussi and Smart

void field_init_curve_ab(field_ptr f, element_ptr a, element_ptr b, mpz_t order, mpz_t cofac);

void field_init_curve_with_map(field_ptr cnew, field_ptr c,
                       field_ptr dstfield, fieldmap map);

void field_init_curve_ab_map(field_t cnew, field_t c,
	fieldmap map, field_ptr mapdest,
	mpz_t ordernew, mpz_t cofacnew);

element_ptr curve_x_coord(element_t e);
element_ptr curve_y_coord(element_t e);
element_ptr curve_a_coeff(element_t e);
element_ptr curve_b_coeff(element_t e);
element_ptr curve_field_a_coeff(field_t f);
element_ptr curve_field_b_coeff(field_t f);

void curve_from_x(element_ptr e, element_t x);
void curve_set_si(element_t R, long int x, long int y);
void curve_set_gen_no_cofac(element_ptr a);

void field_curve_use_random_solvefory(field_ptr f);

#endif //__PBC_CURVE_H__
