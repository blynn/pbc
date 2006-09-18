#include <assert.h>
#include "pbc.h"
#include "param.h"
#include "f_param.h"
#include "fieldquadratic.h"
#include "tracker.h"
#include "utils.h"

struct f_pairing_data_s {
    field_t Fq, Fq2, Fq2x, Fq12;
    curve_t Eq, Etwist;
    element_t alphainv2, alphainv3;
    mpz_t tateexp;
    fieldmap mapbase;
    element_t xpowq2, xpowq6, xpowq8;
};
typedef struct f_pairing_data_s f_pairing_data_t[1];
typedef struct f_pairing_data_s *f_pairing_data_ptr;

void f_param_init(f_param_t fp)
{
    mpz_init(fp->q);
    mpz_init(fp->r);
    mpz_init(fp->b);
    mpz_init(fp->beta);
    mpz_init(fp->alpha0);
    mpz_init(fp->alpha1);
}

void f_param_clear(f_param_t fp)
{
    mpz_clear(fp->q);
    mpz_clear(fp->r);
    mpz_clear(fp->b);
    mpz_clear(fp->beta);
    mpz_clear(fp->alpha0);
    mpz_clear(fp->alpha1);
}

void f_param_out_str(FILE *stream, f_param_ptr p)
{
    param_out_type(stream, "f");
    param_out_mpz(stream, "q", p->q);
    param_out_mpz(stream, "r", p->r);
    param_out_mpz(stream, "b", p->b);
    param_out_mpz(stream, "beta", p->beta);
    param_out_mpz(stream, "alpha0", p->alpha0);
    param_out_mpz(stream, "alpha1", p->alpha1);
}

void f_param_inp_generic (f_param_ptr p, fetch_ops_t *fops, void *ctx)
{
    assert (fops);
    assert (ctx);
    symtab_t tab;

    symtab_init(tab);
    param_read_generic (tab, fops, ctx);

    lookup_mpz(p->q, tab, "q");
    lookup_mpz(p->r, tab, "r");
    lookup_mpz(p->b, tab, "b");
    lookup_mpz(p->beta, tab, "beta");
    lookup_mpz(p->alpha0, tab, "alpha0");
    lookup_mpz(p->alpha1, tab, "alpha1");

    param_clear_tab(tab);
    symtab_clear(tab);
}

void f_param_inp_buf (f_param_ptr p, const char *buf, size_t len)
{
    assert (buf);
    tracker_t t;
    tracker_init (&t, buf, len);
    f_param_inp_generic (p, &fops_buf, &t);
}

void f_param_inp_str (f_param_ptr p, FILE *stream)
{
    assert (stream);
    f_param_inp_generic (p, &fops_str, stream);
}

static void tryminusx(mpz_ptr q, mpz_ptr x)
{
    //36x4 - 36x3 + 24x2 - 6x + 1
    //= ((36(x - 1)x + 24)x - 6)x + 1
    mpz_sub_ui(q, x, 1);
    mpz_mul(q, q, x);
    mpz_mul_ui(q, q, 36);
    mpz_add_ui(q, q, 24);
    mpz_mul(q, q, x);
    mpz_sub_ui(q, q, 6);
    mpz_mul(q, q, x);
    mpz_add_ui(q, q, 1);
}

static void tryplusx(mpz_ptr q, mpz_ptr x)
{
    //36x4 + 36x3 + 24x2 + 6x + 1
    //= ((36(x + 1)x + 24)x + 6)x + 1
    mpz_add_ui(q, x, 1);
    mpz_mul(q, q, x);
    mpz_mul_ui(q, q, 36);
    mpz_add_ui(q, q, 24);
    mpz_mul(q, q, x);
    mpz_add_ui(q, q, 6);
    mpz_mul(q, q, x);
    mpz_add_ui(q, q, 1);
}

void f_param_gen(f_param_t fp, int bits)
{
    //36 is a 6-bit number
    int xbit = (bits - 6) / 4;
    //TODO: use binary search to find smallest appropriate x
    mpz_t x, t;
    mpz_ptr q = fp->q;
    mpz_ptr r = fp->r;
    mpz_ptr b = fp->b;
    field_t Fq, Fq2, Fq2x;
    element_t e1;
    element_t f;
    curve_t c;
    point_t P;

    mpz_init(x);
    mpz_init(t);
    mpz_setbit(x, xbit);
    for (;;) {
	mpz_mul(t, x, x);
	mpz_mul_ui(t, t, 6);
	mpz_add_ui(t, t, 1);
	tryminusx(q, x);
	mpz_sub(r, q, t);
	mpz_add_ui(r, r, 1);
	if (mpz_probab_prime_p(q, 10) && mpz_probab_prime_p(r, 10)) break;

	tryplusx(q, x);
	mpz_sub(r, q, t);
	mpz_add_ui(r, r, 1);
	if (mpz_probab_prime_p(q, 10) && mpz_probab_prime_p(r, 10)) break;

	mpz_add_ui(x, x, 1);
    }

    field_init_fp(Fq, q);
    element_init(e1, Fq);

    for (;;) {
	element_random(e1);
	curve_init_b(c, e1);
	point_init(P, c);

	point_random(P);

	point_mul(P, r, P);
	if (point_is_inf(P)) break;
	point_clear(P);
	curve_clear(c);
    }
    element_to_mpz(b, e1);

    element_clear(e1);

    field_init_quadratic(Fq2, Fq);

    element_to_mpz(fp->beta, field_get_nqr(Fq));

    field_init_poly(Fq2x, Fq2);

    element_init(f, Fq2x);

    //find irreducible polynomial f = x^6 + alpha
    poly_alloc(f, 7);
    element_set1(poly_coeff(f, 6));
    for (;;) {
	element_random(poly_coeff(f, 0));
	if (poly_is_irred(f)) break;
    }

    //extend F_q^2 using f = x^6 + alpha
    //see if sextic twist contains a subgroup of order r
    //if not, it's the wrong twist: replace alpha with alpha^5
    {
	curve_t ctest;
	point_t Ptest;
	mpz_t z0, z1;
	mpz_init(z0);
	mpz_init(z1);
	element_init(e1, Fq2);
	element_set_mpz(e1, fp->b);
	element_mul(e1, e1, poly_coeff(f, 0));
	element_neg(e1, e1);

	curve_init_b(ctest, e1);
	point_init(Ptest, ctest);
	point_random(Ptest);

	//I'm not sure what the #E'(F_q^2) is, but
	//it definitely divides n_12 = #E(F_q^12). It contains a
	//subgroup of order r if and only if
	//(n_12 / r^2)P != O for some (in fact most) P in E'(F_q^6)
	mpz_pow_ui(z0, q, 12);
	mpz_add_ui(z0, z0, 1);
	compute_trace_n(z1, q, t, 12);
	mpz_sub(z1, z0, z1);
	mpz_mul(z0, r, r);
	mpz_divexact(z1, z1, z0);

	point_mul(Ptest, z1, Ptest);
	if (point_is_inf(Ptest)) {
	    mpz_set_ui(z0, 5);
	    element_pow(poly_coeff(f, 0), poly_coeff(f, 0), z0);
	}
	element_clear(e1);
	point_clear(Ptest);
	curve_clear(ctest);
	mpz_clear(z0);
	mpz_clear(z1);
    }

    element_to_mpz(fp->alpha0, fi_re(poly_coeff(f, 0)));
    element_to_mpz(fp->alpha1, fi_im(poly_coeff(f, 0)));

    element_clear(f);

    field_clear(Fq2x);
    field_clear(Fq2);
    field_clear(Fq);

    mpz_clear(t);
    mpz_clear(x);
}

static void Fq_to_Fq12(element_t out, element_t in)
{
    element_set0(out);
    element_set(fi_re(polymod_coeff(out, 0)), in);
}

static void cc_miller_no_denom(element_t res, mpz_t q, point_t P,
	element_ptr Qx, element_ptr Qy, fieldmap mapbase)
{
    int m;
    element_t v;
    point_t Z;
    element_t a, b, c;
    const common_curve_ptr cc = P->curve->data;
    element_t t0;
    element_t e0, e1;
    UNUSED_VAR (cc);

    void do_vertical(void)
    {
	mapbase(e0, Z->x);
	element_sub(e0, Qx, e0);
	element_mul(v, v, e0);
    }

    void do_tangent(void)
    {
	//a = -3 Zx^2 since cc->a is 0 for D = 3
	//b = 2 * Zy
	//c = -(2 Zy^2 + a Zx);
	const element_ptr Zx = Z->x;
	const element_ptr Zy = Z->y;

	element_square(a, Zx);
	element_mul_si(a, a, 3);
	element_neg(a, a);

	element_add(b, Zy, Zy);

	element_mul(t0, b, Zy);
	element_mul(c, a, Zx);
	element_add(c, c, t0);
	element_neg(c, c);

	//TODO: use poly_mul_constant?
	mapbase(e0, a);
	element_mul(e0, e0, Qx);
	mapbase(e1, b);
	element_mul(e1, e1, Qy);
	element_add(e0, e0, e1);
	mapbase(e1, c);
	element_add(e0, e0, e1);
	element_mul(v, v, e0);
    }

    void do_line(void)
    {
	//a = -(B.y - A.y) / (B.x - A.x);
	//b = 1;
	//c = -(A.y + a * A.x);
	//but we'll multiply by B.x - A.x to avoid division

	const element_ptr Ax = Z->x;
	const element_ptr Ay = Z->y;
	const element_ptr Bx = P->x;
	const element_ptr By = P->y;

	element_sub(b, Bx, Ax);
	element_sub(a, Ay, By);
	element_mul(t0, b, Ay);
	element_mul(c, a, Ax);
	element_add(c, c, t0);
	element_neg(c, c);

	mapbase(e0, a);
	element_mul(e0, e0, Qx);
	mapbase(e1, b);
	element_mul(e1, e1, Qy);
	element_add(e0, e0, e1);
	mapbase(e1, c);
	element_add(e0, e0, e1);
	element_mul(v, v, e0);
    }

    element_init(a, P->curve->field);
    element_init(b, P->curve->field);
    element_init(c, P->curve->field);
    element_init(t0, P->curve->field);
    element_init(e0, res->field);
    element_init(e1, res->field);

    element_init(v, res->field);
    point_init(Z, P->curve);

    point_set(Z, P);

    element_set1(v);
    m = mpz_sizeinbase(q, 2) - 2;

    for(;;) {
	do_tangent();

	if (!m) break;

	point_double(Z, Z);
	if (mpz_tstbit(q, m)) {
	    do_line();
	    point_add(Z, Z, P);
	}
	m--;
	element_square(v, v);
    }

    element_set(res, v);

    element_clear(v);
    point_clear(Z);
    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(t0);
    element_clear(e0);
    element_clear(e1);
}

static void f_pairing(element_ptr out, element_ptr in1, element_ptr in2,
	pairing_t pairing)
{
    point_ptr Qbase = in2->data;
    element_t x, y;
    f_pairing_data_ptr p = pairing->data;

    element_init(x, out->field);
    element_init(y, out->field);
    //map from twist: (x, y) --> (v^-2 x, v^-3 y)
    //where v is the sixth root used to construct the twist
    polymod_const_mul(x, Qbase->x, p->alphainv2);
    polymod_const_mul(y, Qbase->y, p->alphainv3);
/*
{
    curve_t c;
    cc_init_map_curve(c, p->Eq, out->field, Fq_to_Fq12);
    point_t P;
    point_init(P, c);
    point_random(P);
    element_set(x, P->x);
    element_set(y, P->y);
}
*/
    cc_miller_no_denom(out, pairing->r, in1->data, x, y, p->mapbase);

{
    element_t epow;
    void qpower(element_ptr e) {
	element_set(polymod_coeff(x, 0), polymod_coeff(out, 0));
	element_mul(polymod_coeff(x, 1), polymod_coeff(out, 1), e);
	element_square(epow, e);
	element_mul(polymod_coeff(x, 2), polymod_coeff(out, 2), epow);
	element_mul(epow, epow, e);
	element_mul(polymod_coeff(x, 3), polymod_coeff(out, 3), epow);
	element_mul(epow, epow, e);
	element_mul(polymod_coeff(x, 4), polymod_coeff(out, 4), epow);
	element_mul(epow, epow, e);
	element_mul(polymod_coeff(x, 5), polymod_coeff(out, 5), epow);
    }
    element_init(epow, p->Fq2);

    qpower(p->xpowq8);
    element_set(y, x);
    qpower(p->xpowq6);
    element_mul(y, y, x);
    qpower(p->xpowq2);
    element_mul(x, x, out);
    element_invert(x, x);
    element_mul(out, y, x);

    element_clear(epow);
}
    element_pow(out, out, p->tateexp);
    element_clear(x);
    element_clear(y);
}

void pairing_init_f_param(pairing_t pairing, f_param_t param)
{
    f_pairing_data_ptr p;
    element_t irred;
    element_t e0, e1, e2;
    mpz_t one;
    p = pairing->data = malloc(sizeof(f_pairing_data_t));
    mpz_init(pairing->r);
    mpz_set(pairing->r, param->r);
    field_init_fp(pairing->Zr, pairing->r);
    field_init_fp(p->Fq, param->q);
    p->Fq->nqr = malloc(sizeof(element_t));
    element_init(p->Fq->nqr, p->Fq);
    element_set_mpz(p->Fq->nqr, param->beta);
    field_init_quadratic(p->Fq2, p->Fq);
    field_init_poly(p->Fq2x, p->Fq2);
    element_init(irred, p->Fq2x);
    poly_alloc(irred, 7);
    element_set1(poly_coeff(irred, 6));
    element_set_mpz(fi_re(poly_coeff(irred, 0)), param->alpha0);
    element_set_mpz(fi_im(poly_coeff(irred, 0)), param->alpha1);
    field_init_polymod(p->Fq12, irred);
    element_clear(irred);

    element_init(e0, p->Fq);
    element_init(e1, p->Fq);
    element_init(e2, p->Fq2);

    element_set_mpz(e1, param->b);
    curve_init_cc_ab(p->Eq, e0, e1);
    element_set_mpz(e0, param->alpha0);
    element_neg(e0, e0);
    element_mul(fi_re(e2), e0, e1);
    element_set_mpz(e0, param->alpha1);
    element_neg(e0, e0);
    element_mul(fi_im(e2), e0, e1);
    element_clear(e0);
    element_init(e0, p->Fq2);
    curve_init_cc_ab(p->Etwist, e0, e2);
    element_clear(e0);
    element_clear(e1);
    element_clear(e2);

    pairing->G1 = malloc(sizeof(field_t));
    pairing->G2 = malloc(sizeof(field_t));

    mpz_init(one);
    mpz_set_si(one, 1);
    field_init_curve_group(pairing->G1, p->Eq, one);
    field_init_curve_group(pairing->G2, p->Etwist, one);
    mpz_clear(one);
    pairing->GT = p->Fq12;
    pairing->map = f_pairing;
    
    element_init(p->alphainv2, p->Fq12);
    element_init(p->alphainv3, p->Fq12);
    element_set1(polymod_coeff(p->alphainv2, 2));
    element_invert(p->alphainv2, p->alphainv2);
    element_set1(polymod_coeff(p->alphainv3, 3));
    element_invert(p->alphainv3, p->alphainv3);
    mpz_init(p->tateexp);
    /*
    mpz_pow_ui(p->tateexp, param->q, 12);
    mpz_sub_ui(p->tateexp, p->tateexp, 1);
    mpz_divexact(p->tateexp, p->tateexp, param->r);
    */
    mpz_ptr z = p->tateexp;
    mpz_mul(z, param->q, param->q);
    mpz_sub_ui(z, z, 1);
    mpz_mul(z, z, param->q);
    mpz_mul(z, z, param->q);
    mpz_add_ui(z, z, 1);
    mpz_divexact(z, z, param->r);

    element_init(p->xpowq2, p->Fq2);
    element_init(p->xpowq6, p->Fq2);
    element_init(p->xpowq8, p->Fq2);
    element_t xpowq;
    element_init(xpowq, p->Fq12);

    //there are smarter ways since we know q = 1 mod 6
    //and that x^6 = -alpha
    //but this is fast enough
    element_set1(polymod_coeff(xpowq, 1));
    element_pow(xpowq, xpowq, param->q);
    element_pow(xpowq, xpowq, param->q);
    element_set(p->xpowq2, polymod_coeff(xpowq, 1));

    element_pow(xpowq, xpowq, param->q);
    element_pow(xpowq, xpowq, param->q);
    element_pow(xpowq, xpowq, param->q);
    element_pow(xpowq, xpowq, param->q);
    element_set(p->xpowq6, polymod_coeff(xpowq, 1));

    element_pow(xpowq, xpowq, param->q);
    element_pow(xpowq, xpowq, param->q);
    element_set(p->xpowq8, polymod_coeff(xpowq, 1));

    element_clear(xpowq);

    p->mapbase = Fq_to_Fq12;
}
