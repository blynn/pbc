// Step-by-step Weil and Tate pairings.
// For my thesis.
#include <string.h>
#include "pbc.h"
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
    /*
    //we could multiply by -2*V.y to avoid division so:
    //a = -(3 Vx^2 + cc->a)
    //b = 2 * Vy
    //c = -(2 Vy^2 + a Vx);
    //
    //actually no, since fasterweil won't work if we do this
    */
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

    element_sub(a, Py, Vy);
    element_sub(b, Vx, Px);
    element_div(a, a, b);
    element_set1(b);
    element_mul(c, a, Vx);
    element_add(c, c, Vy);
    element_neg(c, c);

    /*
    //but we could multiply by B.x - A.x to avoid division, so
    //a = -(By - Ay)
    //b = Bx - Ax
    //c = -(Ay b + a Ax);
    element_sub(a, Vy, Py);
    element_sub(b, Px, Vx);
    element_mul(c, Vx, Py);
    element_mul(e0, Vy, Px);
    element_sub(c, c, e0);
    //
    //actually no, since fasterweil won't work if we do this
    */

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

void millertate(element_t z, element_t P, element_t Q)
{
    element_t Z;
    element_t z0;

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

    element_clear(z0);
    element_clear(Z);
}

void tate(element_t z, element_t P, element_t Q)
{
    mpz_t q1r;

    mpz_init(q1r);
    mpz_set_ui(q1r, 696);

    /*
    millertate(z, P, Q);
    element_printf("prepow: z = %B\n", z);
    element_pow_mpz(z, z, q1r);
    */
    {
        element_t R, QR;
        element_t z0;

        element_init_same_as(R, P);
        element_init_same_as(QR, P);
        element_init_same_as(z0, z);

        element_random(R);
        element_add(QR, Q, R);

        millertate(z, P, QR);
        millertate(z0, P, R);
        element_div(z, z, z0);
        element_pow_mpz(z, z, q1r);
        element_clear(R);
        element_clear(QR);
    }

    mpz_clear(q1r);
}

void shipseystange(element_t z, element_t P, element_t Q)
{
    mpz_t q1r;

    mpz_init(q1r);
    mpz_set_ui(q1r, 696);

    element_ptr x = curve_x_coord(P);
    element_ptr y = curve_y_coord(P);

    element_ptr x2 = curve_x_coord(Q);
    element_ptr y2 = curve_y_coord(Q);

    element_t v0m1, v0m2, v0m3;
    element_t v00, v01, v02, v03, v04;
    element_t v1m1, v10, v11;
    element_t t0, t1, t2;
    element_t W20inv;
    element_t Wm11inv;
    element_t W2m1inv;
    element_t sm2, sm1, s0, s1, s2, s3;
    element_t pm2, pm1, p0, p1, p2, p3;

    element_init_same_as(sm2, z);
    element_init_same_as(sm1, z);
    element_init_same_as(s0, z);
    element_init_same_as(s1, z);
    element_init_same_as(s2, z);
    element_init_same_as(s3, z);

    element_init_same_as(pm2, z);
    element_init_same_as(pm1, z);
    element_init_same_as(p0, z);
    element_init_same_as(p1, z);
    element_init_same_as(p2, z);
    element_init_same_as(p3, z);

    element_init_same_as(v0m3, z);
    element_init_same_as(v0m2, z);
    element_init_same_as(v0m1, z);
    element_init_same_as(v00, z);
    element_init_same_as(v01, z);
    element_init_same_as(v02, z);
    element_init_same_as(v03, z);
    element_init_same_as(v04, z);

    element_init_same_as(v1m1, z);
    element_init_same_as(v10, z);
    element_init_same_as(v11, z);

    element_init_same_as(W20inv, z);
    element_init_same_as(Wm11inv, z);
    element_init_same_as(W2m1inv, z);

    element_init_same_as(t0, z);
    element_init_same_as(t1, z);
    element_init_same_as(t2, z);

    element_set0(v0m1);
    element_set1(v00);
    element_neg(v0m2, v00);
    element_double(v01, y);

    element_neg(v0m3, v01);

    element_invert(W20inv, v01);

    element_sub(Wm11inv, x, x2);
    element_square(t1, Wm11inv);
    element_invert(Wm11inv, Wm11inv);
    element_double(t0, x);
    element_add(t0, t0, x2);
    element_mul(t1, t0, t1);
    element_add(t0, y, y2);
    element_square(t0, t0);
    element_sub(t0, t0, t1);
    element_invert(W2m1inv, t0);

    /* Let P=(x,y) since A=1, B=0 we have:
     * W(3,0) = 3x^4 + 6x^2 - 1
     * W(4,0) = 4y(x^6 + 5x^4 - 5x^2 - 1)
     */

    //t0 = x^2
    element_square(t0, x);

    //t1 = x^4
    element_square(t1, t0);

    //t2 = x^4 + 2 x^2
    element_double(t2, t0);
    element_add(t2, t2, t1);

    //v02 = W(3,0)
    element_double(v02, t2);
    element_add(v02, v02, t2);
    element_add(v02, v02, v0m2);

    //t2 = x^4 - x^2
    element_sub(t2, t1, t0);

    //v03 = 5(x^4 - x^2)
    element_double(v03, t2);
    element_double(v03, v03);
    element_add(v03, v03, t2);

    //t2 = x^6
    element_mul(t2, t0, t1);

    //v03 = W(4,0)
    element_add(v03, v03, t2);
    element_add(v03, v03, v0m2);
    element_double(v03, v03);
    element_double(v03, v03);
    element_mul(v03, v03, y);

    //v04 = W(5,0) = W(2,0)^3 W(4,0) - W(3,0)^3
    element_square(t0, v01);
    element_mul(t0, t0, v01);
    element_mul(v04, t0, v03);
    element_square(t0, v02);
    element_mul(t0, t0, v02);
    element_sub(v04, v04, t0);

    element_set1(v1m1);
    element_set1(v10);

    element_printf("x y: %B %B\n", x, y);
    element_printf("x2 y2: %B %B\n", x2, y2);
    element_sub(t0, x2, x);
    element_sub(t1, y2, y);
    element_div(t0, t1, t0);
    element_square(t0, t0);
    element_double(v11, x);
    element_add(v11, v11, x2);
    element_sub(v11, v11, t0);

    element_printf("VEC1: %B %B %B\n", v1m1, v10, v11);
    element_printf("VEC0: %B %B %B %B %B %B %B %B\n",
            v0m3, v0m2, v0m1, v00, v01, v02, v03, v04);

    //Double
    element_square(sm2, v0m2);
    element_square(sm1, v0m1);
    element_square(s0, v00);
    element_square(s1, v01);
    element_square(s2, v02);
    element_square(s3, v03);

    element_mul(pm2, v0m3, v0m1);
    element_mul(pm1, v0m2, v00);
    element_mul(p0, v0m1, v01);
    element_mul(p1, v00, v02);
    element_mul(p2, v01, v03);
    element_mul(p3, v02, v04);

    element_mul(t0, pm1, sm2);
    element_mul(t1, pm2, sm1);
    element_sub(v0m3, t0, t1);

    element_mul(t1, pm2, s0);
    element_mul(t0, p0, sm2);
    element_sub(v0m2, t0, t1);
    element_mul(v0m2, v0m2, W20inv);

    element_mul(t0, p0, sm1);
    element_mul(t1, pm1, s0);
    element_sub(v0m1, t0, t1);

    element_mul(t1, pm1, s1);
    element_mul(t0, p1, sm1);
    element_sub(v00, t0, t1);
    element_mul(v00, v00, W20inv);

    element_mul(t0, p1, s0);
    element_mul(t1, p0, s1);
    element_sub(v01, t0, t1);

    element_mul(t1, p0, s2);
    element_mul(t0, p2, s0);
    element_sub(v02, t0, t1);
    element_mul(v02, v02, W20inv);

    element_mul(t0, p2, s1);
    element_mul(t1, p1, s2);
    element_sub(v03, t0, t1);

    element_mul(t1, p1, s3);
    element_mul(t0, p3, s1);
    element_sub(v04, t0, t1);
    element_mul(v04, v04, W20inv);

    element_square(t0, v10);
    element_mul(t1, v1m1, v11);

    element_mul(t2, pm1, t0);
    element_mul(v1m1, t1, sm1);
    element_sub(v1m1, v1m1, t2);

    element_mul(t2, p0, t0);
    element_mul(v10, t1, s0);
    element_sub(v10, v10, t2);

    element_mul(t2, p1, t0);
    element_mul(v11, t1, s1);
    element_sub(v11, v11, t2);
    element_mul(v11, v11, Wm11inv);

    element_printf("VEC1: %B %B %B\n", v1m1, v10, v11);
    element_printf("VEC0: %B %B %B %B %B %B %B %B\n",
            v0m3, v0m2, v0m1, v00, v01, v02, v03, v04);

    //DoubleAdd
    element_square(sm2, v0m2);
    element_square(sm1, v0m1);
    element_square(s0, v00);
    element_square(s1, v01);
    element_square(s2, v02);
    element_square(s3, v03);

    element_mul(pm2, v0m3, v0m1);
    element_mul(pm1, v0m2, v00);
    element_mul(p0, v0m1, v01);
    element_mul(p1, v00, v02);
    element_mul(p2, v01, v03);
    element_mul(p3, v02, v04);

    element_mul(t1, pm2, s0);
    element_mul(t0, p0, sm2);
    element_sub(v0m3, t0, t1);
    element_mul(v0m3, v0m3, W20inv);

    element_mul(t0, p0, sm1);
    element_mul(t1, pm1, s0);
    element_sub(v0m2, t0, t1);

    element_mul(t1, pm1, s1);
    element_mul(t0, p1, sm1);
    element_sub(v0m1, t0, t1);
    element_mul(v0m1, v0m1, W20inv);

    element_mul(t0, p1, s0);
    element_mul(t1, p0, s1);
    element_sub(v00, t0, t1);

    element_mul(t1, p0, s2);
    element_mul(t0, p2, s0);
    element_sub(v01, t0, t1);
    element_mul(v01, v01, W20inv);

    element_mul(t0, p2, s1);
    element_mul(t1, p1, s2);
    element_sub(v02, t0, t1);

    element_mul(t1, p1, s3);
    element_mul(t0, p3, s1);
    element_sub(v03, t0, t1);
    element_mul(v03, v03, W20inv);

    element_mul(t0, p3, s2);
    element_mul(t1, p2, s3);
    element_sub(v04, t0, t1);

    element_square(t0, v10);
    element_mul(t1, v1m1, v11);

    element_mul(t2, p0, t0);
    element_mul(v1m1, t1, s0);
    element_sub(v1m1, v1m1, t2);

    element_mul(t2, p1, t0);
    element_mul(v10, t1, s1);
    element_sub(v10, v10, t2);
    element_mul(v10, v10, Wm11inv);

    element_mul(t2, t1, s2);
    element_mul(v11, p2, t0);
    element_sub(v11, v11, t2);
    element_mul(v11, v11, W2m1inv);

    element_printf("VEC1: %B %B %B\n", v1m1, v10, v11);
    element_printf("VEC0: %B %B %B %B %B %B %B %B\n",
            v0m3, v0m2, v0m1, v00, v01, v02, v03, v04);
    element_div(z, v11, v01);
    element_printf("prepow: %B\n", z);

    element_pow_mpz(z, z, q1r);

    mpz_clear(q1r);
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

    element_clear(gr);
    element_clear(r);
    element_clear(hs);
    element_clear(s);
    element_clear(z);
    element_clear(z0);
    element_clear(z1);
}

void fasterweil(element_t w, element_t g, element_t h)
{
    element_t hs;
    element_t s;
    element_t z, z0, z1;

    element_init(z, Fq2);
    element_init(z0, Fq2);
    element_init(z1, Fq2);

    element_init_same_as(hs, h);
    element_init_same_as(s, h);

    element_random(s);
    //point_random always takes the same square root
    //why not take the other one for once?
    element_set_str(s, "[[48,55],[28,51]]", 0);

    element_printf("chose S = %B\n", s);
    element_add(hs, h, s);

    element_printf("Q+S = %B\n", hs);

    millertate(z, g, hs);
    millertate(z0, g, s);
    element_div(z1, z, z0);
    element_printf("num: %B\n", z1);

    miller(w, hs, s, h, g);
    element_printf("denom: %B\n", w);

    element_div(w, z1, w);

    element_clear(z);
    element_clear(z0);
    element_clear(z1);
    element_clear(hs);
    element_clear(s);
}

void fasterweil2(element_t w, element_t g, element_t h)
{
    element_t gr;
    element_t r;
    element_t z, z0, z1;

    element_init(z, Fq2);
    element_init(z0, Fq2);
    element_init(z1, Fq2);

    element_init_same_as(gr, g);
    element_init_same_as(r, g);

    element_random(r);
    //point_random always takes the same square root
    //why not take the other one for once?
    element_set_str(r, "[[48,55],[28,51]]", 0);

    element_printf("chose R = %B\n", r);
    element_add(gr, g, r);

    element_printf("P+R = %B\n", gr);

    miller(w, gr, r, g, h);
    element_printf("num: %B\n", w);

    millertate(z, h, gr);
    millertate(z0, h, r);
    element_div(z1, z, z0);
    element_printf("denom: %B\n", z1);

    element_div(w, w, z1);

    element_clear(z);
    element_clear(z0);
    element_clear(z1);
    element_clear(gr);
    element_clear(r);
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

    fasterweil(w0, g, h);
    element_printf("fasterweil: %B\n", w0);

    element_set1(w1);
    for (i=1; i<6; i++) {
        element_mul(w1, w1, w0);
        element_printf("%d: %B\n", i, w1);
    }

    fasterweil2(w0, g, h);
    element_printf("fasterweil2: %B\n", w0);

    tate(w0, g, h);
    element_printf("tate: %B\n", w0);

    element_set1(w1);
    for (i=1; i<6; i++) {
        element_mul(w1, w1, w0);
        element_printf("%d: %B\n", i, w1);
    }

    shipseystange(w0, g, h);
    element_printf("ss-tate: %B\n", w0);

    element_set1(w1);
    for (i=1; i<6; i++) {
        element_mul(w1, w1, w0);
        element_printf("%d: %B\n", i, w1);
    }
    return 0;
}
