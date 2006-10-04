#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "fops.h"
#include "symtab.h"
#include "darray.h"
#include "field.h"
#include "poly.h"
#include "fp.h"
#include "fieldquadratic.h"
#include "hilbert.h"
#include "mnt.h"
#include "curve.h"
#include "pairing.h"
#include "d_param.h"
#include "param.h"
#include "tracker.h"
#include "utils.h"

struct mnt_pairing_data_s {
    field_t Fq, Fqx, Fqd, Fqk;
    field_t Eq, Etwist;
    element_t nqrinv, nqrinv2;
    mpz_t tateexp;
    int k;
    element_ptr *xpowq;
};
typedef struct mnt_pairing_data_s mnt_pairing_data_t[1];
typedef struct mnt_pairing_data_s *mnt_pairing_data_ptr;

void d_param_init(d_param_ptr param)
{
    mpz_init(param->q);
    mpz_init(param->n);
    mpz_init(param->h);
    mpz_init(param->r);
    mpz_init(param->a);
    mpz_init(param->b);
    mpz_init(param->nk);
    mpz_init(param->hk);
    param->k = 0;
    param->coeff = NULL;
    mpz_init(param->nqr);
}

void d_param_clear(d_param_ptr param)
{
    int d = param->k / 2;
    int i;
    mpz_clear(param->q);
    mpz_clear(param->n);
    mpz_clear(param->h);
    mpz_clear(param->r);
    mpz_clear(param->a);
    mpz_clear(param->b);
    mpz_clear(param->nk);
    mpz_clear(param->hk);
    mpz_clear(param->nqr);
    for (i=0; i<d; i++) {
	mpz_clear(param->coeff[i]);
    }
    free(param->coeff);
}

void d_param_out_str(FILE *stream, d_param_ptr p)
{
    int d = p->k / 2;
    int i;
    char s[80];
    param_out_type(stream, "d");
    param_out_mpz(stream, "q", p->q);
    param_out_mpz(stream, "n", p->n);
    param_out_mpz(stream, "h", p->h);
    param_out_mpz(stream, "r", p->r);
    param_out_mpz(stream, "a", p->a);
    param_out_mpz(stream, "b", p->b);
    param_out_int(stream, "k", p->k);
    param_out_mpz(stream, "nk", p->nk);
    param_out_mpz(stream, "hk", p->hk);
    for (i=0; i<d; i++) {
	sprintf(s, "coeff%d", i);
	param_out_mpz(stream, s, p->coeff[i]);
    }
    param_out_mpz(stream, "nqr", p->nqr);
}

void d_param_inp_generic (d_param_ptr p, fetch_ops_t fops, void *ctx)
{
    assert (fops);
    assert (ctx);
    symtab_t tab;
    char s[80];
    int i, d;

    symtab_init(tab);
    param_read_generic (tab, fops, ctx);

    lookup_mpz(p->q, tab, "q");
    lookup_mpz(p->n, tab, "n");
    lookup_mpz(p->h, tab, "h");
    lookup_mpz(p->r, tab, "r");
    lookup_mpz(p->a, tab, "a");
    lookup_mpz(p->b, tab, "b");
    p->k = lookup_int(tab, "k");
    lookup_mpz(p->nk, tab, "nk");
    lookup_mpz(p->hk, tab, "hk");
    lookup_mpz(p->nqr, tab, "nqr");

    d = p->k / 2;
    p->coeff = realloc(p->coeff, sizeof(mpz_t) * d);
    for (i=0; i<d; i++) {
	sprintf(s, "coeff%d", i);
	mpz_init(p->coeff[i]);
	lookup_mpz(p->coeff[i], tab, s);
    }

    param_clear_tab(tab);
    symtab_clear(tab);
}

static inline void d_miller_evalfn(element_t e0,
	element_t a, element_t b, element_t c,
	element_t Qx, element_t Qy)
{
    //a, b, c are in Fq
    //point Q is (Qx, Qy * sqrt(nqr)) where nqr is used to construct
    //the quadratic field extension Fqk of Fqd
    element_ptr re_out = fi_re(e0);
    element_ptr im_out = fi_im(e0);

    int i;
    int d = polymod_field_degree(re_out->field);
    for (i=0; i<d; i++) {
	element_mul(polymod_coeff(re_out, i), polymod_coeff(Qx, i), a);
	element_mul(polymod_coeff(im_out, i), polymod_coeff(Qy, i), b);
    }
    element_add(polymod_coeff(re_out, 0), polymod_coeff(re_out, 0), c);
}

static void cc_miller_no_denom_proj(element_t res, mpz_t q, element_t P,
	element_ptr Qx, element_ptr Qy)
{
    int m;
    element_t v;
    element_t Z;
    element_t a, b, c;
    element_t t0, t1;
    element_ptr t2 = a, t3 = b, t4 = c;
    element_t e0;
    element_t z, z2;
    element_ptr Zx, Zy;
    const element_ptr curve_a = curve_a_coeff(P);
    const element_ptr Px = curve_x_coord(P);
    const element_ptr Py = curve_y_coord(P);

    void proj_double(void)
    {
	//t0 = 3x^2 + (curve_a) z^4
	element_square(t0, Zx);
	//element_mul_si(t0, t0, 3);
	element_double(t1, t0);
	element_add(t0, t0, t1);
	element_square(t1, z2);
	element_mul(t1, t1, curve_a);
	element_add(t0, t0, t1);

	//z_out = 2 y z
	element_mul(z, Zy, z);
	//element_mul_si(z, z, 2);
	element_double(z, z);
	element_square(z2, z);

	//t1 = 4 x y^2
	element_square(t2, Zy);
	element_mul(t1, Zx, t2);
	//element_mul_si(t1, t1, 4);
	element_double(t1, t1);
	element_double(t1, t1);

	//x_out = t0^2 - 2 t1
	//element_mul_si(t3, t1, 2);
	element_double(t3, t1);
	element_square(Zx, t0);
	element_sub(Zx, Zx, t3);

	//t2 = 8y^4
	element_square(t2, t2);
	//element_mul_si(t2, t2, 8);
	element_double(t2, t2);
	element_double(t2, t2);
	element_double(t2, t2);

	//y_out = t0(t1 - x_out) - t2
	element_sub(t1, t1, Zx);
	element_mul(t0, t0, t1);
	element_sub(Zy, t0, t2);
    }

    void proj_mixin(void)
    {
	//t2 = Px z^2
	element_mul(t2, z2, Px);

	//t3 = Zx - t2
	element_sub(t3, Zx, t2);

	//t0 = Py z^3
	element_mul(t0, z2, Py);
	element_mul(t0, t0, z);

	//t1 = Zy - t0
	element_sub(t1, Zy, t0);

	//e7 = Zx + t2, use t2 to double for e7
	element_add(t2, Zx, t2);

	//e8 = Zy + t0, use t0 to double for e8
	element_add(t0, Zy, t0);

	//z = z t3
	element_mul(z, z, t3);
	element_square(z2, z);

	//Zx = t1^2 - e7 t3^2
	//t3 now holds t3^3,
	//t4 holds e7 t3^2
	element_square(t4, t3);
	element_mul(t3, t4, t3);
	element_square(Zx, t1);
	element_mul(t4, t2, t4);
	element_sub(Zx, Zx, t4);

	//t4 = e7 t3^2 - 2 Zx
	element_sub(t4, t4, Zx);
	element_sub(t4, t4, Zx);

	//Zy = (t4 t1 - e8 t3^3)/2
	element_mul(t4, t4, t1);
	element_mul(t0, t0, t3);
	element_sub(t4, t4, t0);
	element_halve(Zy, t4);
    }

    void do_tangent(void)
    {
	//a = -(3x^2 + cca z^4)
	//b = 2 y z^3
	//c = -(2 y^2 + x a)
	//a = z^2 a
	element_square(a, z2);
	element_mul(a, a, curve_a);
	element_square(b, Zx);
	//element_mul_si(b, b, 3);
	element_double(t0, b);
	element_add(b, b, t0);
	element_add(a, a, b);
	element_neg(a, a);

	element_mul(b, z, z2);
	element_mul(b, b, Zy);
	element_mul_si(b, b, 2);

	element_mul(c, Zx, a);
	element_mul(a, a, z2);
	element_square(t0, Zy);
	element_mul_si(t0, t0, 2);
	element_add(c, c, t0);
	element_neg(c, c);

	d_miller_evalfn(e0, a, b, c, Qx, Qy);
	element_mul(v, v, e0);
    }

    void do_line(void)
    {
	//a = -(Py z^3 - Zy)
	//b = Px z^3 - Zx z
	//c = Zx z Py - Zy Px;

	element_mul(t0, Zx, z);
	element_mul(t1, z2, z);

	element_mul(a, Py, t1);
	element_sub(a, Zy, a);

	element_mul(b, Px, t1);
	element_sub(b, b, t0);

	element_mul(t0, t0, Py);
	element_mul(c, Zy, Px);
	element_sub(c, t0, c);

	d_miller_evalfn(e0, a, b, c, Qx, Qy);
	element_mul(v, v, e0);
    }

    element_init(a, Px->field);
    element_init(b, a->field);
    element_init(c, a->field);
    element_init(t0, a->field);
    element_init(t1, a->field);
    element_init(e0, res->field);
    element_init(z, a->field);
    element_init(z2, a->field);
    element_set1(z);
    element_set1(z2);

    element_init(v, res->field);
    element_init(Z, P->field);

    element_set(Z, P);
    Zx = curve_x_coord(Z);
    Zy = curve_x_coord(Z);

    element_set1(v);
    m = mpz_sizeinbase(q, 2) - 2;

    for(;;) {
	do_tangent();
	if (!m) break;
	proj_double();
	if (mpz_tstbit(q, m)) {
	    do_line();
	    proj_mixin();
	}
	m--;
	element_square(v, v);
    }

    element_set(res, v);

    element_clear(v);
    element_clear(Z);
    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(t0);
    element_clear(t1);
    element_clear(e0);
    element_clear(z);
    element_clear(z2);
}

static void cc_miller_no_denom_affine(element_t res, mpz_t q, element_t P,
	element_ptr Qx, element_ptr Qy)
{
    int m;
    element_t v;
    element_t Z;
    element_t a, b, c;
    element_t t0;
    element_t e0;
    const element_ptr cca = curve_a_coeff(P);
    const element_ptr Px = curve_x_coord(P);
    const element_ptr Py = curve_y_coord(P);
    element_ptr Zx, Zy;

    /* TODO: when exactly is this not needed?
    void do_vertical(void)
    {
	mapbase(e0, Z->x);
	element_sub(e0, Qx, e0);
	element_mul(v, v, e0);
    }
    */

    void do_tangent(void)
    {
	//a = -(3 Zx^2 + cc->a)
	//b = 2 * Zy
	//c = -(2 Zy^2 + a Zx);

	element_square(a, Zx);
	element_mul_si(a, a, 3);
	element_add(a, a, cca);
	element_neg(a, a);

	element_add(b, Zy, Zy);

	element_mul(t0, b, Zy);
	element_mul(c, a, Zx);
	element_add(c, c, t0);
	element_neg(c, c);

	d_miller_evalfn(e0, a, b, c, Qx, Qy);
	element_mul(v, v, e0);
    }

    void do_line(void)
    {
	//a = -(B.y - A.y) / (B.x - A.x);
	//b = 1;
	//c = -(A.y + a * A.x);
	//but we'll multiply by B.x - A.x to avoid division

	element_sub(b, Px, Zx);
	element_sub(a, Zy, Py);
	element_mul(t0, b, Zy);
	element_mul(c, a, Zx);
	element_add(c, c, t0);
	element_neg(c, c);

	d_miller_evalfn(e0, a, b, c, Qx, Qy);
	element_mul(v, v, e0);
    }

    element_init(a, Px->field);
    element_init(b, a->field);
    element_init(c, a->field);
    element_init(t0, a->field);
    element_init(e0, res->field);

    element_init(v, res->field);
    element_init(Z, P->field);

    element_set(Z, P);
    Zx = curve_x_coord(Z);
    Zy = curve_y_coord(Z);

    element_set1(v);
    m = mpz_sizeinbase(q, 2) - 2;

    for(;;) {
	do_tangent();

	if (!m) break;

	element_double(Z, Z);
	if (mpz_tstbit(q, m)) {
	    do_line();
	    element_add(Z, Z, P);
	}
	m--;
	element_square(v, v);
    }

    element_set(res, v);

    element_clear(v);
    element_clear(Z);
    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(t0);
    element_clear(e0);
}

static void cc_tatepower(element_ptr out, element_ptr in, pairing_t pairing)
{
    mnt_pairing_data_ptr p = pairing->data;
    if (p->k == 6) {
	element_t e0, e1, e2, e3;
	element_init(e0, p->Fqk);
	element_init(e1, p->Fqd);
	element_init(e2, p->Fqd);
	element_init(e3, p->Fqk);
	element_ptr e0re = fi_re(e0);
	element_ptr e0im = fi_im(e0);
	element_t *inre = fi_re(in)->data;
	element_t *inim = fi_im(in)->data;
	//if beta is quadratic nonresidue in F_q^d, note that
	//beta^q = -beta
	//hence for x^q, x^{q^3} we have to negate the beta part of e0
	//as the following routine does not do this
	void qpower(element_ptr e) {
	    element_set0(e0);
	    element_set(((element_t *) e0re->data)[0], inre[0]);
	    element_set(((element_t *) e0im->data)[0], inim[0]);
	    polymod_const_mul(e2, inre[1], e);
	    element_add(e0re, e0re, e2);
	    polymod_const_mul(e2, inim[1], e);
	    element_add(e0im, e0im, e2);
	    element_square(e1, e);
	    polymod_const_mul(e2, inre[2], e1);
	    element_add(e0re, e0re, e2);
	    polymod_const_mul(e2, inim[2], e1);
	    element_add(e0im, e0im, e2);
	}
	qpower(p->xpowq[4]);
	element_set(e3, e0);
	qpower(p->xpowq[3]);
	element_neg(e0im, e0im);
	element_mul(e3, e3, e0);
	qpower(p->xpowq[1]);
	element_neg(e0im, e0im);
	element_mul(e0, e0, in);
	element_invert(e0, e0);
	element_mul(out, e3, e0);
	element_pow_mpz(out, out, p->tateexp);
	element_clear(e0);
	element_clear(e1);
	element_clear(e2);
	element_clear(e3);
    } else {
	element_pow_mpz(out, in, p->tateexp);
    }
}

static void (*cc_miller_no_denom_fn)(element_t res, mpz_t q, element_t P,
	element_ptr Qx, element_ptr Qy);

static void cc_pairing(element_ptr out, element_ptr in1, element_ptr in2,
	pairing_t pairing)
{
    element_ptr Qbase = in2;
    element_t Qx, Qy;
    mnt_pairing_data_ptr p = pairing->data;

    element_init(Qx, p->Fqd);
    element_init(Qy, p->Fqd);
    //map from twist: (x, y) --> (v^-1 x, v^-(3/2) y)
    //where v is the quadratic nonresidue used to construct the twist
    element_mul(Qx, curve_x_coord(Qbase), p->nqrinv);
    //v^-3/2 = v^-2 * v^1/2
    element_mul(Qy, curve_y_coord(Qbase), p->nqrinv2);
    cc_miller_no_denom_fn(out, pairing->r, in1, Qx, Qy);
    cc_tatepower(out, out, pairing);
    element_clear(Qx);
    element_clear(Qy);
}

static int cc_is_almost_coddh(element_ptr a, element_ptr b,
	element_ptr c, element_ptr d,
	pairing_t pairing)
{
    int res = 0;
    element_t t0, t1, t2;
    element_t cx, cy;
    element_t dx, dy;
    mnt_pairing_data_ptr p = pairing->data;

    element_init(cx, p->Fqd);
    element_init(cy, p->Fqd);
    element_init(dx, p->Fqd);
    element_init(dy, p->Fqd);

    element_init(t0, pairing->GT);
    element_init(t1, pairing->GT);
    element_init(t2, pairing->GT);
    //map from twist: (x, y) --> (v^-1 x, v^-(3/2) y)
    //where v is the quadratic nonresidue used to construct the twist
    element_mul(cx, curve_x_coord(c), p->nqrinv);
    element_mul(dx, curve_x_coord(d), p->nqrinv);
    //v^-3/2 = v^-2 * v^1/2
    element_mul(cy, curve_y_coord(c), p->nqrinv2);
    element_mul(dy, curve_y_coord(d), p->nqrinv2);

    cc_miller_no_denom_fn(t0, pairing->r, a->data, dx, dy);
    cc_miller_no_denom_fn(t1, pairing->r, b->data, cx, cy);
    cc_tatepower(t0, t0, pairing);
    cc_tatepower(t1, t1, pairing);
    element_mul(t2, t0, t1);
    if (element_is1(t2)) {
	//g, g^x, h, h^-x case
	res = 1;
    } else {
	element_invert(t1, t1);
	element_mul(t2, t0, t1);
	if (element_is1(t2)) {
	    //g, g^x, h, h^x case
	    res = 1;
	}
    }
    element_clear(cx);
    element_clear(cy);
    element_clear(dx);
    element_clear(dy);
    element_clear(t0);
    element_clear(t1);
    element_clear(t2);
    return res;
}

static void trace(element_t out, element_t in, pairing_ptr pairing)
{
    UNUSED_VAR(out); UNUSED_VAR(in); UNUSED_VAR(pairing);
    /* TODO: broken
    int i;
    element_ptr p = in;
    element_t r;
    element_t q;
    element_ptr outp = out;
    mnt_pairing_data_ptr pdp = pairing->data;

    element_init(q, pdp->Eqk);
    element_init(r, pdp->Eqk);

    q->inf_flag = 0;
    //map from twist: (x, y) --> (v^-1 x, v^-(3/2) y)
    //where v is the quadratic nonresidue used to construct the twist
    element_mul(fi_re(q->x), p->x, pdp->nqrinv);
    //v^-3/2 = v^-2 * v^1/2
    element_mul(fi_im(q->y), p->y, pdp->nqrinv2);
    point_set(r, q);

    printf("r is ");
    point_out_str(stdout, 0, r);
    printf("\n");

    for (i=1; i<pdp->k; i++) {
	cc_frobenius(q, q, pdp->Fq->order);
    printf("q is ");
    point_out_str(stdout, 0, q);
    printf("\n");

	point_add(r, r, q);
    printf("r is ");
    point_out_str(stdout, 0, r);
    printf("\n");
    }
    point_clear(q);

    printf("r is \n");
    point_out_str(stdout, 0, r);
    printf("\n");
    exit(1);
    element_set(outp->x, polymod_coeff(fi_re(r->x), 0));
    element_set(outp->y, polymod_coeff(fi_re(r->y), 0));
    //may have to multiply to ensure order is correct
    point_clear(r);
    */
}

static void d_pairing_option_set(pairing_t pairing, char *key, char *value)
{
    UNUSED_VAR(pairing);
    if (!strcmp(key, "coord")) {
	if (!strcmp(value, "projective")) {
	    cc_miller_no_denom_fn = cc_miller_no_denom_proj;
	} else if (!strcmp(value, "affine")) {
	    cc_miller_no_denom_fn = cc_miller_no_denom_affine;
	}
    }
}

struct pp_coeff_s {
    element_t a;
    element_t b;
    element_t c;
};
typedef struct pp_coeff_s pp_coeff_t[1];
typedef struct pp_coeff_s *pp_coeff_ptr;

static void d_pairing_pp_init(pairing_pp_t p, element_ptr in1, pairing_t pairing)
{
    element_ptr P = in1;
    const element_ptr Px = curve_x_coord(P);
    const element_ptr Py = curve_y_coord(P);
    element_t Z;
    int m;
    mnt_pairing_data_ptr info = pairing->data;
    element_t t0;
    element_t a, b, c;
    field_ptr Fq = info->Fq;
    pp_coeff_t *coeff;
    mpz_ptr q = pairing->r;
    pp_coeff_ptr pp;
    const element_ptr cca = curve_a_coeff(P);
    element_ptr Zx;
    element_ptr Zy;

    void store_abc(void)
    {
	element_init(pp->a, Fq);
	element_init(pp->b, Fq);
	element_init(pp->c, Fq);
	element_set(pp->a, a);
	element_set(pp->b, b);
	element_set(pp->c, c);
	pp++;
    }

    void do_tangent(void)
    {
	//a = -slope_tangent(Z.x, Z.y);
	//b = 1;
	//c = -(Z.y + a * Z.x);
	//but we multiply by 2*Z.y to avoid division

	//a = -Zx * (3 Zx + twicea_2) - a_4;
	//Common curves: a2 = 0 (and cc->a is a_4), so
	//a = -(3 Zx^2 + cc->a)
	//b = 2 * Zy
	//c = -(2 Zy^2 + a Zx);

	element_square(a, Zx);
	element_double(t0, a);
	element_add(a, a, t0);
	element_add(a, a, cca);
	element_neg(a, a);

	element_add(b, Zy, Zy);

	element_mul(t0, b, Zy);
	element_mul(c, a, Zx);
	element_add(c, c, t0);
	element_neg(c, c);

	store_abc();
    }

    void do_line(void)
    {
	//a = -(B.y - A.y) / (B.x - A.x);
	//b = 1;
	//c = -(A.y + a * A.x);
	//but we'll multiply by B.x - A.x to avoid division

	element_sub(b, Px, Zx);
	element_sub(a, Zy, Py);
	element_mul(t0, b, Zy);
	element_mul(c, a, Zx);
	element_add(c, c, t0);
	element_neg(c, c);

	store_abc();
    }

    element_init(Z, P->field);
    element_set(Z, P);
    Zx = curve_x_coord(Z);
    Zy = curve_y_coord(Z);

    element_init(t0, Fq);
    element_init(a, Fq);
    element_init(b, Fq);
    element_init(c, Fq);

    m = mpz_sizeinbase(q, 2) - 2;
    p->data = malloc(sizeof(pp_coeff_t) * 2 * m);
    coeff = (pp_coeff_t *) p->data;
    pp = coeff[0];

    for(;;) {
	do_tangent();

	if (!m) break;

	element_double(Z, Z);
	if (mpz_tstbit(q, m)) {
	    do_line();
	    element_add(Z, Z, P);
	}
	m--;
    }

    element_clear(t0);
    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(Z);
}

static void d_pairing_pp_clear(pairing_pp_t p)
{
    //TODO: better to store a sentinel value in p->data?
    mpz_ptr q = p->pairing->r;
    int m = mpz_sizeinbase(q, 2) + mpz_popcount(q) - 3;
    int i;
    pp_coeff_t *coeff = (pp_coeff_t *) p->data;
    pp_coeff_ptr pp;
    for (i=0; i<m; i++) {
	pp = coeff[i];
	element_clear(pp->a);
	element_clear(pp->b);
	element_clear(pp->c);
    }
    free(p->data);
}

static void d_pairing_pp_apply(element_ptr out, element_ptr in2, pairing_pp_t p)
{
    mpz_ptr q = p->pairing->r;
    mnt_pairing_data_ptr info = p->pairing->data;
    int m = mpz_sizeinbase(q, 2) - 2;
    pp_coeff_t *coeff = (pp_coeff_t *) p->data;
    pp_coeff_ptr pp = coeff[0];
    element_ptr Qbase = in2;
    element_t e0;
    element_t Qx, Qy;
    element_t v;
    element_init_GT(e0, p->pairing);
    element_init_GT(v, p->pairing);
    element_init(Qx, info->Fqd);
    element_init(Qy, info->Fqd);

    //map from twist: (x, y) --> (v^-1 x, v^-(3/2) y)
    //where v is the quadratic nonresidue used to construct the twist
    element_mul(Qx, curve_x_coord(Qbase), info->nqrinv);
    //v^-3/2 = v^-2 * v^1/2
    element_mul(Qy, curve_y_coord(Qbase), info->nqrinv2);

    element_set1(out);
    for(;;) {
	d_miller_evalfn(e0, pp->a, pp->b, pp->c, Qx, Qy);
	element_mul(out, out, e0);
	pp++;

	if (!m) break;

	if (mpz_tstbit(q, m)) {
	    d_miller_evalfn(e0, pp->a, pp->b, pp->c, Qx, Qy);
	    element_mul(out, out, e0);
	    pp++;
	}
	m--;
	element_square(out, out);
    }
    cc_tatepower(out, out, p->pairing);

    element_clear(e0);
    element_clear(Qx);
    element_clear(Qy);
    element_clear(v);
}

void d_pairing_clear(pairing_t pairing)
{
    mnt_pairing_data_ptr p = pairing->data;

    mpz_clear(p->tateexp);
    if (p->k == 6) {
	int i;
	for (i=1; i<=4; i++) element_clear(p->xpowq[i]);
	free(p->xpowq);
    }

    field_clear(p->Etwist);
    field_clear(p->Eq);
    element_clear(p->nqrinv);
    element_clear(p->nqrinv2);
    field_clear(p->Fqk);
    field_clear(p->Fqd);
    field_clear(p->Fqx);
    field_clear(p->Fq);
}

void pairing_init_d_param(pairing_t pairing, d_param_t param)
{
    mnt_pairing_data_ptr p;
    element_t a, b;
    element_t irred;
    int d = param->k / 2;
    int i;

    if (param->k % 2) {
	fprintf(stderr, "odd k not implemented anymore\n");
	exit(2);
    }

    mpz_init(pairing->r);
    mpz_set(pairing->r, param->r);
    field_init_fp(pairing->Zr, pairing->r);
    pairing->map = cc_pairing;
    pairing->is_almost_coddh = cc_is_almost_coddh;

    p =	pairing->data = malloc(sizeof(mnt_pairing_data_t));
    field_init_fp(p->Fq, param->q);
    element_init(a, p->Fq);
    element_init(b, p->Fq);
    element_set_mpz(a, param->a);
    element_set_mpz(b, param->b);
    field_init_curve_ab(p->Eq, a, b, pairing->r, param->h);

    field_init_poly(p->Fqx, p->Fq);
    element_init(irred, p->Fqx);
    poly_alloc(irred, d + 1);
    for (i=0; i<d; i++) {
	element_set_mpz(poly_coeff(irred, i), param->coeff[i]);
    }
    element_set1(poly_coeff(irred, d));

    field_init_polymod(p->Fqd, irred);
    element_clear(irred);

    p->Fqd->nqr = malloc(sizeof(element_t));
    element_init(p->Fqd->nqr, p->Fqd);
    element_set_mpz(((element_t *) p->Fqd->nqr->data)[0], param->nqr);

    field_init_quadratic(p->Fqk, p->Fqd);

    if (param->k == 6) {
	element_t e0;
	mpz_ptr q = param->q;
	mpz_ptr z = p->tateexp;
	mpz_init(z);
	mpz_mul(z, q, q);
	mpz_sub(z, z, q);
	mpz_add_ui(z, z, 1);
	mpz_divexact(z, z, pairing->r);

	p->xpowq = malloc(sizeof(element_ptr) * 5);
	element_init(e0, p->Fqd);
	element_set1(((element_t *) e0->data)[1]);
	for (i=1; i<=4; i++) {
	    element_ptr e = p->xpowq[i] = malloc(sizeof(element_t));
	    element_init(e, p->Fqd);
	    element_pow_mpz(e0, e0, q);
	    element_set(e, e0);
	}
	element_clear(e0);
    } else {
	mpz_init(p->tateexp);
	mpz_sub_ui(p->tateexp, p->Fqk->order, 1);
	mpz_divexact(p->tateexp, p->tateexp, pairing->r);
    }

    field_init_curve_ab_map(p->Etwist, p->Eq, element_field_to_polymod, p->Fqd, pairing->r, NULL);

    twist_curve(p->Etwist);

    element_init(p->nqrinv, p->Fqd);
    element_invert(p->nqrinv, field_get_nqr(p->Fqd));
    element_init(p->nqrinv2, p->Fqd);
    element_square(p->nqrinv2, p->nqrinv);

    pairing->G1 = p->Eq;
    pairing->G2 = p->Etwist;

    p->k = param->k;
    pairing->GT = p->Fqk;
    pairing->phi = trace;

    cc_miller_no_denom_fn = cc_miller_no_denom_affine;
    pairing->option_set = d_pairing_option_set;
    pairing->pp_init = d_pairing_pp_init;
    pairing->pp_clear = d_pairing_pp_clear;
    pairing->pp_apply = d_pairing_pp_apply;

    pairing->clear_func = d_pairing_clear;

    element_clear(a);
    element_clear(b);
}

static void compute_cm_curve(d_param_ptr param, cm_info_ptr cm)
    //computes a curve and sets fp to the field it is defined over
    //using the complex multiplication method, where cm holds
    //the appropriate information (e.g. discriminant, field order)
{
    darray_t coefflist;
    element_t hp, root;
    field_t fp, fpx;
    int i, n;
    field_t cc;

    field_init_fp(fp, cm->q);
    field_init_poly(fpx, fp);
    element_init(hp, fpx);

    darray_init(coefflist);

    hilbert_poly(coefflist, cm->D);

    n = coefflist->count;
    poly_alloc(hp, n);
    for (i=0; i<n; i++) {
	element_set_mpz(poly_coeff(hp, i), coefflist->item[i]);
    }

    hilbert_poly_clear(coefflist);

    darray_clear(coefflist);
    //TODO: remove x = 0, 1728 roots
    //TODO: what if there's no roots?
    //printf("hp ");
    //element_out_str(stdout, 0, hp);
    //printf("\n");

    element_init(root, fp);
    findroot(root, hp);
    //printf("root = ");
    //element_out_str(stdout, 0, root);
    //printf("\n");
    element_clear(hp);
    field_clear(fpx);

    //the root is the j-invariant of our desired curve
    field_init_curve_j(cc, root, cm->n, NULL);
    element_clear(root);

    //we may need to twist it however
    {
	element_t P;

	//pick a random point P and see if it has the right order
	element_init(P, cc);
	element_random(P);
	element_mul_mpz(P, P, cm->n);
	//element_printf("P = %B", P);
	//if not, we twist the curve
	if (!element_is0(P)) {
	    twist_curve(cc);
	}
	element_clear(P);
    }

    mpz_set(param->q, cm->q);
    mpz_set(param->n, cm->n);
    mpz_set(param->h, cm->h);
    mpz_set(param->r, cm->r);
    element_to_mpz(param->a, curve_field_a_coeff(cc));
    element_to_mpz(param->b, curve_field_b_coeff(cc));
    param->k = cm->k;
    {
	mpz_t z;
	mpz_init(z);
	//compute order of curve in F_q^k
	//n = q - t + 1 hence t = q - n + 1
	mpz_sub(z, param->q, param->n);
	mpz_add_ui(z, z, 1);
	compute_trace_n(z, param->q, z, param->k);
	mpz_pow_ui(param->nk, param->q, param->k);
	mpz_sub_ui(z, z, 1);
	mpz_sub(param->nk, param->nk, z);
	mpz_mul(z, param->r, param->r);
	mpz_divexact(param->hk, param->nk, z);
	mpz_clear(z);
    }
    field_clear(cc);
    field_clear(fp);
}

void d_param_from_cm(d_param_t param, cm_info_ptr cm)
{
    field_t Fq, Fqx, Fqd;
    element_t irred, nqr;
    int d = cm->k / 2;
    int i;

    compute_cm_curve(param, cm);

    field_init_fp(Fq, param->q);
    field_init_poly(Fqx, Fq);
    element_init(irred, Fqx);
    do {
	poly_random_monic(irred, d);
    } while (!poly_is_irred(irred));
    field_init_polymod(Fqd, irred);

    //find a quadratic nonresidue of Fqd lying in Fq
    element_init(nqr, Fqd);
    do {
	element_random(((element_t *) nqr->data)[0]);
    } while (element_is_sqr(nqr));

    param->coeff = realloc(param->coeff, sizeof(mpz_t) * d);

    for (i=0; i<d; i++) {
	mpz_init(param->coeff[i]);
	element_to_mpz(param->coeff[i], poly_coeff(irred, i));
    }
    element_to_mpz(param->nqr, ((element_t *) nqr->data)[0]);

    element_clear(nqr);
    element_clear(irred);

    field_clear(Fqx);
    field_clear(Fqd);
    field_clear(Fq);
}
