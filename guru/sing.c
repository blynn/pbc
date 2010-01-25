/*
 * Example of a singular curve, similar to 19.c
 * but the Tate pairing degenerates
 *
 * Consider the curve E: y^2 = x^3 + x^2 over F_19:
 * E_ns(F_19) is a cyclic group of order 18.
 */

#include "pbc.h"
#include "pbc_singular.h"
#include "pbc_fp.h"

static void miller(element_t res, element_t P, element_t Q, element_t R, int n)
{
    //collate divisions
    int m;
    element_t v, vd;
    element_t Z;
    element_t a, b, c;
    element_t e0, e1;
    mpz_t q;
    element_ptr Zx, Zy;
    const element_ptr Px = curve_x_coord(P);
    const element_ptr Py = curve_y_coord(P);
    const element_ptr numx = curve_x_coord(Q);
    const element_ptr numy = curve_y_coord(Q);
    const element_ptr denomx = curve_x_coord(R);
    const element_ptr denomy = curve_y_coord(R);

    void do_vertical(element_t e, element_t edenom)
    {
        element_sub(e0, numx, Zx);
        element_mul(e, e, e0);

        element_sub(e0, denomx, Zx);
        element_mul(edenom, edenom, e0);
    }

    void do_tangent(element_t e, element_t edenom)
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

        if (element_is0(Zy)) {
            do_vertical(e, edenom);
            return;
        }
        element_square(a, Zx);
        element_mul_si(a, a, 3);
        element_add(a, a, Zx);
        element_add(a, a, Zx);
        element_neg(a, a);

        element_add(b, Zy, Zy);

        element_mul(e0, b, Zy);
        element_mul(c, a, Zx);
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

    void do_line(element_ptr e, element_ptr edenom)
    {
        if (!element_cmp(Zx, Px)) {
            if (!element_cmp(Zy, Py)) {
                do_tangent(e, edenom);
            } else {
                do_vertical(e, edenom);
            }
            return;
        }

        element_sub(b, Px, Zx);
        element_sub(a, Zy, Py);
        element_mul(c, Zx, Py);
        element_mul(e0, Zy, Px);
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
    element_init(Z, P->field);

    element_set(Z, P);
    Zx = curve_x_coord(Z);
    Zy = curve_y_coord(Z);

    element_set1(v);
    element_set1(vd);

    mpz_init(q);
    mpz_set_ui(q, n);
    m = mpz_sizeinbase(q, 2) - 2;

    while(m >= 0) {
        element_square(v, v);
        element_square(vd, vd);
        do_tangent(v, vd);
        element_double(Z, Z);
        do_vertical(vd, v);

        if (mpz_tstbit(q, m)) {
            do_line(v, vd);
            element_add(Z, Z, P);
            if (m) {
                do_vertical(vd, v);
            }
        }
        m--;
    }

    mpz_clear(q);

    element_invert(vd, vd);
    element_mul(res, v, vd);

    element_clear(v);
    element_clear(vd);
    element_clear(Z);
    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(e0);
    element_clear(e1);
}

static void tate_3(element_ptr out, element_ptr P, element_ptr Q, element_ptr R)
{
    mpz_t six;

    mpz_init(six);
    mpz_set_ui(six, 6);
    element_t QR;
    element_t e0;

    element_init(QR, P->field);
    element_init(e0, out->field);

    element_add(QR, Q, R);

    //for subgroup size 3, -2P = P, hence
    //the tangent line at P has divisor 3(P) - 3(O)

    miller(out, P, QR, R, 3);

    element_pow_mpz(out, out, six);
    element_clear(QR);
    element_clear(e0);
    mpz_clear(six);
}

static void tate_9(element_ptr out, element_ptr P, element_ptr Q, element_ptr R)
{
    element_t QR;
    element_init(QR, P->field);

    element_add(QR, Q, R);

    miller(out, P, QR, R, 9);

    element_square(out, out);

    element_clear(QR);
}

int main(void)
{
    field_t c;
    field_t Z19;
    element_t P, Q, R;
    mpz_t q, z;
    element_t a;
    int i;

    mpz_init(q);
    mpz_init(z);

    mpz_set_ui(q, 19);

    field_init_fp(Z19, q);
    element_init(a, Z19);

    field_init_curve_singular_with_node(c, Z19);

    element_init(P, c);
    element_init(Q, c);
    element_init(R, c);

    //(3,+/-6) is a generator
    //we have an isomorphism from E_ns to F_19^*
    // (3,6) --> 3
    //(generally (x,y) --> (y+x)/(y-x)

    curve_set_si(R, 3, 6);

    for (i=1; i<=18; i++) {
        mpz_set_si(z, i);
        element_mul_mpz(Q, R, z);
        element_printf("%dR = %B\n", i, Q);
    }

    mpz_set_ui(z, 6);
    element_mul_mpz(P, R, z);
    //P has order 3
    element_printf("P = %B\n", P);

    for (i=1; i<=3; i++) {
        mpz_set_si(z, i);
        element_mul_mpz(Q, R, z);
        tate_3(a, P, Q, R);
        element_printf("e_3(P,%dP) = %B\n", i, a);
    }

    element_double(P, R);
    //P has order 9
    element_printf("P = %B\n", P);
    for (i=1; i<=9; i++) {
        mpz_set_si(z, i);
        element_mul_mpz(Q, P, z);
        tate_9(a, P, Q, R);
        element_printf("e_9(P,%dP) = %B\n", i, a);
    }

    return 0;
}
