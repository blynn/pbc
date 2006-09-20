/* 
 * Example of a singular curve, similar to 19.c
 * but the Tate pairing degenerates
 *
 * Consider the curve E: y^2 = x^3 + x^2 over F_19:
 * E_ns(F_19) is a cyclic group of order 18.
 */

#include "pbc.h"

static void miller(element_t res, point_t P,
	element_ptr numx, element_ptr numy,
	element_ptr denomx, element_ptr denomy, int n)
{
    //collate divisions
    int m;
    element_t v, vd;
    point_t Z;
    element_t a, b, c;
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
	//This curve is special:
	//a = -(3 Ax^2 + 2Ax)
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
	element_add(a, a, Ax);
	element_add(a, a, Ax);
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

int main(void)
{
    curve_t c;
    field_t Z19;
    point_t P, Q, R;
    mpz_t q, z;
    element_t a;
    int i;

    mpz_init(q);
    mpz_init(z);

    mpz_set_ui(q, 19);

    field_init_fp(Z19, q);
    element_init(a, Z19);

    curve_init_singular_with_node(c, Z19);

    point_init(P, c);
    point_init(Q, c);
    point_init(R, c);

    //(3,+/-6) is a generator
    //we have an isomorphism from E_ns to F_19^*
    // (3,6) --> 3
    //(generally (x,y) --> (y+x)/(y-x)

    //TODO: write wrappers for this sort of thing
    element_set_si(R->x, 3);
    element_set_si(R->y, 6);
    R->inf_flag = 0;

    for (i=1; i<=18; i++) {
	mpz_set_si(z, i);
	point_mul(Q, z, R);
	printf("%dR = ", i);
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
	printf("e_3(P,%dP) = ", i);
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
	point_mul(Q, z, P);
	tate_9(a, P, Q, R);
	printf("e_9(P,%dP) = ", i);
	element_out_str(stdout, 0, a);
	printf("\n");
    }

    return 0;
}
