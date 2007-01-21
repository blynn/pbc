//step-by-step Weil and Tate pairings
//for my thesis
//Ben Lynn
#include <string.h>
#include <pbc.h>
#include "pbc_fp.h"
#include "pbc_fieldquadratic.h"

static field_t Fq, Fq2, E, E2;
static mpz_t order;

static void do_vert(element_ptr z, element_ptr V, element_ptr Q)
{
    element_ptr Vx = curve_x_coord(V);
    element_ptr Qx = curve_x_coord(Q);
    element_ptr Qy = curve_y_coord(Q);

    element_t a, b, c;
    element_init_same_as(a, Vx);
    element_init_same_as(b, Vx);
    element_init_same_as(c, Vx);

    //a = 1
    //b = 0;
    //c = -Vx
    element_set1(a);
    element_set0(b);
    element_neg(c, Vx);

    element_printf("vert at %B: %B %B %B\n", Vx, a, b, c);
    element_mul(a, a, Qx);
    element_mul(b, b, Qy);
    element_add(c, c, a);
    element_add(z, c, b);
    element_printf("vert eval = %B\n", z);
    element_clear(a);
    element_clear(b);
    element_clear(c);
}

static void do_tangent(element_ptr z, element_ptr V, element_ptr Q)
{
    element_ptr Vx = curve_x_coord(V);
    element_ptr Vy = curve_y_coord(V);
    element_ptr Qx = curve_x_coord(Q);
    element_ptr Qy = curve_y_coord(Q);

    element_t a, b, c;
    element_init_same_as(a, Vx);
    element_init_same_as(b, Vx);
    element_init_same_as(c, Vx);

    //a = -slope_tangent(V.x, V.y);
    //b = 1;
    //c = -(V.y + aV.x);
    //we could multiply by -2*V.y to avoid division so:
    //a = -(3 Vx^2 + cc->a)
    //b = 2 * Vy
    //c = -(2 Vy^2 + a Vx);
    element_square(a, Vx);
    //element_mul_si(a, a, 3);
    element_add(b, a, a);
    element_add(a, b, a);
    element_set1(b);
    element_add(a, a, b);
    element_neg(a, a);
    element_double(b, Vy);
    element_div(a, a, b);
    element_set1(b);
    element_mul(c, a, Vx);
    element_add(c, c, Vy);
    element_neg(c, c);

    element_printf("tan at %B: %B %B %B\n", V, a, b, c);

    element_mul(a, a, Qx);
    element_mul(b, b, Qy);
    element_add(c, c, a);
    element_add(z, c, b);
    element_printf("tan eval = %B\n", z);
    element_clear(a);
    element_clear(b);
    element_clear(c);
}

static void do_line(element_ptr z, element_ptr V, element_ptr P, element_ptr Q)
{
    element_ptr Vx = curve_x_coord(V);
    element_ptr Vy = curve_y_coord(V);
    element_ptr Px = curve_x_coord(P);
    element_ptr Py = curve_y_coord(P);
    element_ptr Qx = curve_x_coord(Q);
    element_ptr Qy = curve_y_coord(Q);

    element_t a, b, c, e0;
    element_init_same_as(a, Vx);
    element_init_same_as(b, Vx);
    element_init_same_as(c, Vx);
    element_init_same_as(e0, Vx);

    //a = -(B.y - A.y) / (B.x - A.x);
    //b = 1;
    //c = -(A.y + a * A.x);
    //but we'll multiply by B.x - A.x to avoid division, so
    //a = -(By - Ay)
    //b = Bx - Ax
    //c = -(Ay b + a Ax);
    element_sub(a, Vy, Py);
    element_sub(b, Px, Vx);
    element_mul(c, Vx, Py);
    element_mul(e0, Vy, Px);
    element_sub(c, c, e0);

    element_printf("line at %B: %B %B %B\n", V, a, b, c);
    element_mul(a, a, Qx);
    element_mul(b, b, Qy);
    element_add(c, c, a);
    element_add(z, c, b);
    element_printf(" = %B\n", z);

    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(e0);
}

void tate(element_t z, element_t P, element_t Q)
{
    element_t Z;
    element_t z0;
    mpz_t q1r;

    mpz_init(q1r);
    element_init_same_as(Z, P);
    element_init_same_as(z0, z);

    element_set(Z, P);

    do_tangent(z, Z, Q);

    element_double(Z, Z);

    do_vert(z0, Z, Q);
    element_div(z, z, z0);

    element_printf("presquare: z = %B\n", z);

    element_square(z, z);

    element_printf("square: z = %B\n", z);

    do_tangent(z0, Z, Q);
    element_mul(z, z, z0);

    mpz_set_ui(q1r, 696);

    element_printf("prepow: z = %B\n", z);
    element_pow_mpz(z, z, q1r);

    mpz_clear(q1r);
    element_clear(z0);
    element_clear(Z);
}

void miller(element_t z, element_t PR, element_t R, element_t P, element_t Q)
{
    int m = mpz_sizeinbase(order, 2) - 2;

    element_t Z;
    element_t z1;
    element_t x1;
    element_init_same_as(Z, PR);

    element_set(Z, P);
    element_set1(z);
    element_init_same_as(z1, z);
    element_init_same_as(x1, z);

    do_vert(x1, PR, Q);
    element_printf("vert(P+R) %B\n", x1);
    do_line(z1, P, R, Q);
    element_printf("line(P,R) %B\n", z1);
    element_div(x1, x1, z1);
    element_printf("x1 %B\n", x1);
    element_set(z, x1);

    for (;;) {
	printf("iteration %d: %d\n", m, mpz_tstbit(order,m));
	element_square(z, z);
	element_printf("squared: %B\n", z);
	do_tangent(z1, Z, Q);
	element_mul(z, z, z1);

	element_double(Z, Z);
	do_vert(z1, Z, Q);
	element_div(z, z, z1);
	element_printf("pre-if: %B\n", z);

	if (mpz_tstbit(order, m)) {
	    element_mul(z, z, x1);
	    do_vert(z1, P, Q);
	    element_mul(z, z, z1);
	    element_printf("done %B\n", z);
	    /*
	    do_line(z1, Z, P, Q);
	    element_mul(z, z, z1);
	    element_add(Z, Z, P);
	    do_vert(z1, Z, Q);
	    element_div(z, z, z1);
	    */
	}
	if (!m) break;
	m--;
    }

    element_clear(x1);
    element_clear(z1);
}
/*
*/

void weil(element_t w, element_t g, element_t h)
{
    element_t gr;
    element_t hs;
    element_t r;
    element_t s;
    element_t z, z0, z1;

    element_init(z, Fq2);
    element_init(z0, Fq2);
    element_init(z1, Fq2);

    element_init_same_as(gr, g);
    element_init_same_as(hs, h);
    element_init_same_as(r, g);
    element_init_same_as(s, h);

    element_random(r);
    element_random(s);
    //point_random always takes the same square root
    //why not take the other one for once?
    element_neg(r, r);
    element_set_str(r, "[[40,0],[54,0]]", 0);
    element_set_str(s, "[[48,55],[28,51]]", 0);

    element_printf("chose R = %B\n", r);
    element_printf("chose S = %B\n", s);
    element_add(gr, g, r);
    element_add(hs, h, s);

    element_printf("P+R = %B\n", gr);
    element_printf("Q+S = %B\n", hs);
    miller(z, gr, r, g, hs);
    miller(z0, gr, r, g, s);
    element_div(z1, z, z0);
    element_printf("num: %B\n", z1);

    miller(z, hs, s, h, gr);
    miller(z0, hs, s, h, r);
    element_div(w, z, z0);
    element_printf("denom: %B\n", w);

    element_div(w, z1, w);
}

int main(void)
{
    int i;
    element_t g, h;
    element_t w0, w1;
    element_t a, b;
    mpz_t prime, cofac;

    mpz_init(prime);
    mpz_init(order);
    mpz_init(cofac);
    mpz_set_ui(prime, 59);

    field_init_fp(Fq, prime);

    element_init(a, Fq);
    element_init(b, Fq);

    field_init_fi(Fq2, Fq);

    element_set1(a);
    element_set0(b);
    mpz_set_ui(order, 5);
    mpz_set_ui(cofac, 12);

    field_init_curve_ab(E, a, b, order, cofac);

    element_clear(a);
    element_clear(b);
    element_init(a, Fq2);
    element_init(b, Fq2);
    element_set1(a);
    element_set0(b);

    mpz_mul(cofac, cofac, cofac);
    field_init_curve_ab(E2, a, b, order, NULL);

    element_init(g, E2);
    element_init(h, E2);

    element_init(w0, Fq2);
    element_init(w1, Fq2);

    /*
    do {
	element_random(g);
    } while (element_is1(g));
    for (i=1; i<5; i++) {
	element_mul(h, h, g);
	element_printf("%d: %B\n", i, h);
	element_printf("tangent = ");
	do_tangent(h);
    }
    */
    element_set_str(g, "[[25,0],[30,0]", 0);
    element_set_str(h, "[[34,0],[0,30]", 0);
    weil(w0, g, h);
    element_printf("weil: %B\n", w0);

    element_set1(w1);
    for (i=1; i<6; i++) {
	element_mul(w1, w1, w0);
	element_printf("%d: %B\n", i, w1);
    }

    tate(w0, g, h);
    element_printf("tate: %B\n", w0);

    element_set1(w1);
    for (i=1; i<6; i++) {
	element_mul(w1, w1, w0);
	element_printf("%d: %B\n", i, w1);
    }
    return 0;
}
