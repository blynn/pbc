#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "symtab.h"
#include "fops.h"
#include "field.h"
#include "fp.h"
#include "fieldquadratic.h"
#include "pairing.h"
#include "a1_param.h"
#include "curve.h"
#include "param.h"
#include "tracker.h"

struct a1_pairing_data_s {
    field_t Fp, Fp2, Ep;
    mpz_t h;
};
typedef struct a1_pairing_data_s a1_pairing_data_t[1];
typedef struct a1_pairing_data_s *a1_pairing_data_ptr;

void a1_param_init(a1_param_t param)
{
    mpz_init(param->p);
    mpz_init(param->n);
}

void a1_param_clear(a1_param_t param)
{
    mpz_clear(param->p);
    mpz_clear(param->n);
}

void a1_param_out_str(FILE *stream, a1_param_ptr p)
{
    param_out_type(stream, "a1");
    param_out_mpz(stream, "p", p->p);
    param_out_mpz(stream, "n", p->n);
    param_out_int(stream, "l", p->l);
}

void a1_param_inp_generic (a1_param_ptr p, fetch_ops_t fops, void *ctx)
{
    assert (fops);
    assert (ctx);
    symtab_t tab;

    symtab_init(tab);
    param_read_generic (tab, fops, ctx);

    lookup_mpz(p->p, tab, "p");
    lookup_mpz(p->n, tab, "n");
    p->l = lookup_int(tab, "l");

    param_clear_tab(tab);
    symtab_clear(tab);
}

void a1_param_gen(a1_param_t param, mpz_t order)
{
    //if order is even, ideally check all even l, not just multiples of 4
    //but I don't see a good reason for having an even order
    unsigned int l = 4;
    mpz_t n;
    mpz_ptr p = param->p;
    mpz_init(n);
    mpz_mul_ui(n, order, 4);
    mpz_sub_ui(p, n, 1);
    for (;;) {
	if (mpz_probab_prime_p(p, 20)) {
	    break;
	}
	mpz_add(p, p, n);
	l += 4;
    }
    param->l = l;
    mpz_set(param->n, order);
    mpz_clear(n);
}

static void phi_identity(element_ptr out, element_ptr in, pairing_ptr pairing)
{
    (void) pairing;
    element_set(out, in);
}

static void a1_pairing(element_ptr out, element_ptr in1, element_ptr in2,
	pairing_t pairing)
//in1, in2 are from E(F_q), out from F_q^2
{
    a1_pairing_data_ptr p = pairing->data;
    element_t V;
    element_t f, f0;
    element_t a, b, c;
    element_t e0;
    int m;
    element_ptr Qx = curve_x_coord(in2);
    element_ptr Qy = curve_y_coord(in2);
    element_ptr Vx;
    element_ptr Vy;

    void do_tangent(void) {
	//a = -slope_tangent(V.x, V.y);
	//b = 1;
	//c = -(V.y + aV.x);
	//but we multiply by -2*V.y to avoid division so:
	//a = -(3 Vx^2 + cc->a)
	//b = 2 * Vy
	//c = -(2 Vy^2 + a Vx);
	element_square(a, Vx);
	element_mul_si(a, a, 3);
	element_set1(b);
	element_add(a, a, b);
	element_neg(a, a);

	element_add(b, Vy, Vy);

	element_mul(e0, b, Vy);
	element_mul(c, a, Vx);
	element_add(c, c, e0);
	element_neg(c, c);

	//we'll map Q via (x,y) --> (-x, iy)
	//hence a Qx + c = -a Qx + c is real while
	//(b Qy) = b Qy i is purely imaginary.
	element_mul(a, a, Qx);
	element_sub(fi_re(f0), c, a);
	element_mul(fi_im(f0), b, Qy);
	element_mul(f, f, f0);
    }

    //TODO: this only called once, get rid of arguments
    void do_line(element_ptr A, element_ptr B) {
	//a = -(B.y - A.y) / (B.x - A.x);
	//b = 1;
	//c = -(A.y + a * A.x);
	//but we'll multiply by B.x - A.x to avoid division, so
	//a = -(By - Ay)
	//b = Bx - Ax
	//c = -(Ay b + a Ax);
	element_ptr Ax = curve_x_coord(A);
	element_ptr Ay = curve_y_coord(A);
	element_ptr Bx = curve_x_coord(B);
	element_ptr By = curve_y_coord(B);

	element_sub(a, Ay, By);
	element_sub(b, Bx, Ax);
	element_mul(c, Ax, By);
	element_mul(e0, Ay, Bx);
	element_sub(c, c, e0);

	//we'll map Q via (x,y) --> (-x, iy)
	//hence a Qx + c = -a Qx + c is real while
	//(b Qy) = b Qy i is purely imaginary.
	element_mul(a, a, Qx);
	element_sub(fi_re(f0), c, a);
	element_mul(fi_im(f0), b, Qy);
	element_mul(f, f, f0);
    }

    element_init(V, p->Ep);
    element_set(V, in1);
    Vx = curve_x_coord(V);
    Vy = curve_y_coord(V);

    element_init(f, p->Fp2);
    element_init(f0, p->Fp2);
    element_set1(f);
    element_init(a, p->Fp);
    element_init(b, p->Fp);
    element_init(c, p->Fp);
    element_init(e0, p->Fp);

    m = mpz_sizeinbase(pairing->r, 2) - 2;

    //sliding NAF
    while(m >= 0) {
	do_tangent();
	if (!m) break;

	element_double(V, V);
	if (mpz_tstbit(pairing->r, m)) {
	    do_line(V, in1);
	    element_add(V, V, in1);
	}

	m--;
	element_square(f, f);
    }

    //Tate exponentiation
    //simpler but slower:
    //element_pow_mpz(out, f, p->tateexp);
    //use this trick instead:
    element_invert(f0, f);
    element_neg(fi_im(f), fi_im(f));
    element_mul(f, f, f0);
    element_pow_mpz(out, f, p->h);

    element_clear(f);
    element_clear(f0);
    element_clear(V);
    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(e0);
}

void a1_pairing_clear(pairing_t pairing)
{
    a1_pairing_data_ptr p = pairing->data;
    mpz_clear(p->h);
    field_clear(p->Ep);
    field_clear(p->Fp2);
    field_clear(p->Fp);
    free(p);

    mpz_clear(pairing->r);
    field_clear(pairing->Zr);
}

void pairing_init_a1_param(pairing_t pairing, a1_param_t param)
{
    element_t a, b;
    mpz_init(pairing->r);
    mpz_set(pairing->r, param->n);
    field_init_fp(pairing->Zr, pairing->r);

    a1_pairing_data_ptr p;

    p =	pairing->data = malloc(sizeof(a1_pairing_data_t));
    mpz_init(p->h);
    mpz_set_ui(p->h, param->l);

    field_init_fp(p->Fp, param->p);
    element_init(a, p->Fp);
    element_init(b, p->Fp);
    element_set1(a);
    element_set0(b);
    field_init_curve_ab(p->Ep, a, b, pairing->r, p->h);
    element_clear(a);
    element_clear(b);
    field_init_fi(p->Fp2, p->Fp);

    pairing->G1 = malloc(sizeof(field_t));
    pairing->G2 = pairing->G1 = p->Ep;
    //pairing->phi = phi_identity;
    pairing->GT = p->Fp2;

    pairing->map = a1_pairing;
    pairing->phi = phi_identity;

    pairing->clear_func = a1_pairing_clear;
}
