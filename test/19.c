/* 
 * Toy example of a field where the Tate pairing can be used
 * but the Weil pairing cannot.
 *
 * Consider the curve E: y^2 = x^3 + x + 6 over F_19:
 * E(F_19) is a cyclic group of order 18.
 * Thus E[3] is not contained in F_19
 * (it turns out E[3] is contained in F_19^3).
 *
 * Hence the Weil pairing cannot be defined over F_19
 * However, F_19 contains the cube roots of unity
 * so we can compute the Tate pairing
 */

/*
 * P = (12,13) generates a group of order 3:
 * <(12,13)> = {(12,13), (12,6), O}
 * e(P,P) = 7, so we have the isomorphism
 * <(12,13)> = <7> (in F_19^*)
 *
 * Similarly P = (4, 6) generates a group of order 9, and we find
 * <(4,6)> = <4>
 *
 * P = (0, 5) generates all of E(F_19)
 * Miller's algorithm will not allow us to calculate e(P, P) without
 * first extending F_19.
 * Instead of extending, we could manipulate rational functions since
 * 19 is small enough that an explicit expression of f_P can be found.
 */

#include "pbc.h"
#include "fp.h"
#include "fieldquadratic.h"

static void miller(element_t res, point_t P,
	element_ptr numx, element_ptr numy,
	element_ptr denomx, element_ptr denomy, int n)
{
    //collate divisions
    int m;
    element_t v, vd;
    point_t Z;
    element_t a, b, c;
    const common_curve_ptr cc = P->curve->data;
    element_t e0, e1;
    mpz_t q;

    void do_vertical(element_t e, element_t edenom, point_ptr A)
    {
	element_sub(e0, numx, A->x);
	element_mul(e, e, e0);

	element_sub(e0, denomx, A->x);
	element_mul(edenom, edenom, e0);
    }

    void do_tangent(element_t e, element_t edenom, point_ptr A)
    {
	//a = -slope_tangent(A.x, A.y);
	//b = 1;
	//c = -(A.y + a * A.x);
	//but we multiply by 2*A.y to avoid division

	//a = -Ax * (Ax + Ax + Ax + twicea_2) - a_4;
	//Common curves: a2 = 0 (and cc->a is a_4), so
	//a = -(3 Ax^2 + cc->a)
	//b = 2 * Ay
	//c = -(2 Ay^2 + a Ax);
	const element_ptr Ax = A->x;
	const element_ptr Ay = A->y;

	if (element_is0(Ay)) {
	    do_vertical(e, edenom, A);
	    return;
	}
	element_square(a, Ax);
	element_mul_si(a, a, 3);
	element_add(a, a, cc->a);
	element_neg(a, a);

	element_add(b, Ay, Ay);

	element_mul(e0, b, Ay);
	element_mul(c, a, Ax);
	element_add(c, c, e0);
	element_neg(c, c);

	element_mul(e0, a, numx);
	element_mul(e1, b, numy);
	element_add(e0, e0, e1);
	element_add(e0, e0, c);
	element_mul(e, e, e0);

	element_mul(e0, a, denomx);
	element_mul(e1, b, denomy);
	element_add(e0, e0, e1);
	element_add(e0, e0, c);
	element_mul(edenom, edenom, e0);
    }

    void do_line(element_ptr e, element_ptr edenom, point_ptr A, point_ptr B)
    {
	const element_ptr Ax = A->x;
	const element_ptr Ay = A->y;
	const element_ptr Bx = B->x;
	const element_ptr By = B->y;

	if (!element_cmp(Ax, Bx)) {
	    if (!element_cmp(Ay, By)) {
		do_tangent(e, edenom, A);
	    } else {
		do_vertical(e, edenom, A);
	    }
	    return;
	}

	element_sub(b, Bx, Ax);
	element_sub(a, Ay, By);
	element_mul(c, Ax, By);
	element_mul(e0, Ay, Bx);
	element_sub(c, c, e0);

	element_mul(e0, a, numx);
	element_mul(e1, b, numy);
	element_add(e0, e0, e1);
	element_add(e0, e0, c);
	element_mul(e, e, e0);

	element_mul(e0, a, denomx);
	element_mul(e1, b, denomy);
	element_add(e0, e0, e1);
	element_add(e0, e0, c);
	element_mul(edenom, edenom, e0);
    }

    element_init(a, res->field);
    element_init(b, res->field);
    element_init(c, res->field);
    element_init(e0, res->field);
    element_init(e1, res->field);

    element_init(v, res->field);
    element_init(vd, res->field);
    point_init(Z, P->curve);

    point_set(Z, P);

    element_set1(v);
    element_set1(vd);

    mpz_init(q);
    mpz_set_ui(q, n);
    m = mpz_sizeinbase(q, 2) - 2;

    while(m >= 0) {
	element_square(v, v);
	element_square(vd, vd);
	do_tangent(v, vd, Z);
	point_double(Z, Z);
	do_vertical(vd, v, Z);

	if (mpz_tstbit(q, m)) {
	    do_line(v, vd, Z, P);
	    point_add(Z, Z, P);
	    if (m) {
		do_vertical(vd, v, Z);
	    }
	}
	m--;
    }

    mpz_clear(q);

    element_invert(vd, vd);
    element_mul(res, v, vd);

    element_clear(v);
    element_clear(vd);
    point_clear(Z);
    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(e0);
    element_clear(e1);
}

static void tate_3(element_ptr out, point_ptr P, point_ptr Q, point_ptr R)
{
    mpz_t six;

    mpz_init(six);
    mpz_set_ui(six, 6);
    point_t QR;
    element_t e0;

    point_init(QR, P->curve);
    element_init(e0, out->field);

    point_add(QR, Q, R);

    //for subgroup size 3, -2P = P, hence
    //the tangent line at P has divisor 3(P) - 3(O)

    miller(out, P, QR->x, QR->y, R->x, R->y, 3);

    element_pow_mpz(out, out, six);
    point_clear(QR);
    element_clear(e0);
    mpz_clear(six);
}

static void tate_9(element_ptr out, point_ptr P, point_ptr Q, point_ptr R)
{
    point_t QR;
    point_init(QR, P->curve);

    point_add(QR, Q, R);

    miller(out, P, QR->x, QR->y, R->x, R->y, 9);

    element_square(out, out);

    point_clear(QR);
}

static void tate_18(element_ptr out, point_ptr P, point_ptr Q, point_ptr R, point_ptr S)
{
    mpz_t pow;
    point_t PR;
    point_t QS;
    point_init(PR, P->curve);
    point_init(QS, P->curve);
    element_t outd;

    element_init(outd, out->field);

    mpz_init(pow);
    mpz_set_ui(pow, (19*19-1)/18);

    point_add(PR, P, R);
    point_add(QS, Q, S);

    if (point_is_inf(QS)) {
	point_t S2;
	point_init(S2, P->curve);
	point_double(S2, S);
	miller(out, PR, S->x, S->y, S2->x, S2->y, 18);
	miller(outd, R, S->x, S->y, S2->x, S2->y, 18);
	point_clear(S2);
    } else {
	miller(out, PR, QS->x, QS->y, S->x, S->y, 18);
	miller(outd, R, QS->x, QS->y, S->x, S->y, 18);
    }

    point_clear(PR);
    point_clear(QS);

    element_invert(outd, outd);
    element_mul(out, out, outd);
    element_pow_mpz(out, out, pow);

    element_clear(outd);
    mpz_clear(pow);
}

int main(void)
{
    curve_t c;
    field_t Z19;
    point_t P, Q, R;
    mpz_t q, z;
    element_t a, b;
    int i;

    field_t Z19_2;
    curve_t c2;
    point_t P2, Q2, R2;
    element_t a2;

    mpz_init(q);
    mpz_init(z);

    mpz_set_ui(q, 19);

    field_init_fp(Z19, q);
    element_init(a, Z19);
    element_init(b, Z19);

    element_set_si(a, 1);
    element_set_si(b, 6);

    curve_init_cc_ab(c, a, b);
    point_init(P, c);
    point_init(Q, c);
    point_init(R, c);

    printf("Y^2 = X^3 + X + 6 over F_19\n");
    //(0,+/-5) is a generator
    element_set0(a);
    point_from_x(R, a);

    for (i=1; i<19; i++) {
	printf("%dR = ", i);
	mpz_set_si(z, i);
	point_mul(Q, z, R);
	point_out_str(stdout, 0, Q);
	printf("\n");
    }

    mpz_set_ui(z, 6);
    point_mul(P, z, R);
    //P has order 3
    printf("P = ");
    point_out_str(stdout, 0, P);
    printf("\n");

    for (i=1; i<=3; i++) {
	mpz_set_si(z, i);
	point_mul(Q, z, R);
	tate_3(a, P, Q, R);
	printf("e_3(P,%dR) = ", i);
	element_out_str(stdout, 0, a);
	printf("\n");
    }

    point_double(P, R);
    //P has order 9
    printf("P = ");
    point_out_str(stdout, 0, P);
    printf("\n");
    for (i=1; i<=9; i++) {
	mpz_set_si(z, i);
	//we're supposed to use multiples of R
	//but 2R works just as well and it allows us
	//to use R as the offset every time
	point_mul(Q, z, P);
	tate_9(a, P, Q, R);
	printf("e_9(P,%dP) = ", i);
	element_out_str(stdout, 0, a);
	printf("\n");
    }

    //to do the pairing on all of E(F_19) we need to move to F_19^2
    //or compute the rational function explicitly
    printf("moving to F_19^2\n");
    field_init_fi(Z19_2, Z19);

    cc_init_map_curve(c2, c, Z19_2, element_field_to_fi);
    point_init(P2, c2);
    point_init(Q2, c2);
    point_init(R2, c2);

    element_init(a2, Z19_2);
    element_set0(a2);
    point_from_x(P2, a2);

    point_random(R2);

    printf("P = ");
    point_out_str(stdout, 0, P2);
    printf("\n");

    for (i=1; i<=18; i++) {
	mpz_set_si(z, i);
	point_mul(Q2, z, P2);
	tate_18(a2, P2, Q2, R2, P2);
	printf("e_18(P,%dP) = ", i);
	element_out_str(stdout, 0, a2);
	printf("\n");
    }

    point_clear(P2);
    point_clear(Q2);
    point_clear(R2);
    element_clear(a2);
    curve_clear(c2);
    field_clear(Z19_2);

    curve_clear(c);
    element_clear(a);
    element_clear(b);
    point_clear(P);
    point_clear(Q);
    point_clear(R);
    field_clear(Z19);

    mpz_clear(q);
    mpz_clear(z);
    return 0;
}
