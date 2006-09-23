#include <assert.h>
#include "e_param.h"
#include "param.h"
#include "tracker.h"

struct e_pairing_data_s {
    field_t Fq;
    curve_t Eq;
    mpz_t tateexp;
    int exp2, exp1;
    int sign1, sign0;
};
typedef struct e_pairing_data_s e_pairing_data_t[1];
typedef struct e_pairing_data_s *e_pairing_data_ptr;

void e_param_init(e_param_t ep)
{
    mpz_init(ep->q);
    mpz_init(ep->r);
    mpz_init(ep->h);
    mpz_init(ep->a);
    mpz_init(ep->b);
}

void e_param_clear(e_param_t ep)
{
    mpz_clear(ep->q);
    mpz_clear(ep->r);
    mpz_clear(ep->h);
    mpz_clear(ep->a);
    mpz_clear(ep->b);
}

void e_param_gen(e_param_t p, int rbits, int qbits)
{
    //3 takes 2 bits to represent
    int hbits = (qbits - 2) / 2 - rbits;
    mpz_ptr q = p->q;
    mpz_ptr r = p->r;
    mpz_ptr h = p->h;
    mpz_t n;
    field_t Fq;
    curve_t cc;
    element_t j;
    int found = 0;

    mpz_init(n);

    do {
	int i;
	mpz_set_ui(r, 0);

	if (rand() % 2) {
	    p->exp2 = rbits - 1;
	    p->sign1 = 1;
	} else {
	    p->exp2 = rbits;
	    p->sign1 = -1;
	}
	mpz_setbit(r, p->exp2);

	p->exp1 = (rand() % (p->exp2 - 1)) + 1;
	//use q as a temp variable
	mpz_set_ui(q, 0);
	mpz_setbit(q, p->exp1);

	if (p->sign1 > 0) {
	    mpz_add(r, r, q);
	} else {
	    mpz_sub(r, r, q);
	}

	if (rand() % 2) {
	    p->sign0 = 1;
	    mpz_add_ui(r, r, 1);
	} else {
	    p->sign0 = -1;
	    mpz_sub_ui(r, r, 1);
	}
	if (!mpz_probab_prime_p(r, 10)) continue;
	for (i=0; i<10; i++) {
	    //use q as a temp variable
	    mpz_set_ui(q, 0);
	    mpz_setbit(q, hbits + 1);
	    pbc_mpz_random(h, q);
	    mpz_mul(h, h, h);
	    mpz_mul_ui(h, h, 3);
	    //finally q takes the value it should
	    mpz_mul(n, r, r);
	    mpz_mul(n, n, h);
	    mpz_add_ui(q, n, 1);
	    if (mpz_probab_prime_p(q, 10)) {
		found = 1;
		break;
	    }
	}
    } while (!found);
    /*
    do {
	mpz_set_ui(r, 0);
	mpz_setbit(r, rbits);
	pbc_mpz_random(r, r);
	mpz_nextprime(r, r);
	mpz_mul(n, r, r);
	mpz_mul_ui(n, n, 3);
	mpz_add_ui(q, n, 1);
    } while (!mpz_probab_prime_p(q, 10));
    */

    field_init_fp(Fq, q);
    /*
    element_init(j, Fq);
    element_set_si(j, -3375);
    curve_init_cc_j(cc, j);
    element_clear(j);
    */
    element_init(j, Fq);
    element_set_si(j, 1);
    curve_init_b(cc, j);
    element_clear(j);

    //we may need to twist it however
    {
	point_t P;

	//pick a random point P and see if it has the right order
	point_init(P, cc);
	point_random(P);
	point_mul(P, n, P);
	//if not, we twist the curve
	if (!point_is_inf(P)) {
	    twist_curve(cc);
	}
	point_clear(P);
    }
    mpz_set(p->a, ((common_curve_ptr) cc->data)->a->data);
    mpz_set(p->b, ((common_curve_ptr) cc->data)->b->data);

    mpz_clear(n);
}

void e_param_out_str(FILE *stream, e_param_ptr p)
{
    param_out_type(stream, "e");
    param_out_mpz(stream, "q", p->q);
    param_out_mpz(stream, "r", p->r);
    param_out_mpz(stream, "h", p->h);
    param_out_mpz(stream, "a", p->a);
    param_out_mpz(stream, "b", p->b);
    param_out_int(stream, "exp2", p->exp2);
    param_out_int(stream, "exp1", p->exp1);
    param_out_int(stream, "sign1", p->sign1);
    param_out_int(stream, "sign0", p->sign0);
}

void e_param_inp_generic (e_param_ptr p, fetch_ops_t fops, void *ctx)
{
    assert (fops);
    assert (ctx);
    symtab_t tab;

    symtab_init(tab);
    param_read_generic (tab, fops, ctx);

    lookup_mpz(p->q, tab, "q");
    lookup_mpz(p->r, tab, "r");
    lookup_mpz(p->h, tab, "h");
    lookup_mpz(p->a, tab, "a");
    lookup_mpz(p->b, tab, "b");
    p->exp2 = lookup_int(tab, "exp2");
    p->exp1 = lookup_int(tab, "exp1");
    p->sign1 = lookup_int(tab, "sign1");
    p->sign0 = lookup_int(tab, "sign0");

    param_clear_tab(tab);
    symtab_clear(tab);
}

static void e_miller_proj(element_t res, point_t P,
	element_ptr numx, element_ptr numy,
	element_ptr denomx, element_ptr denomy,
	e_pairing_data_ptr p)
{
    //collate divisions
    int n;
    element_t v, vd;
    element_t v1, vd1;
    point_t Z, Z1;
    element_t a, b, c;
    const common_curve_ptr cc = P->curve->data;
    const element_ptr cca = cc->a;
    element_t e0, e1;
    const element_ptr e2 = a, e3 = b;
    element_t z, z2;
    int i;

    //convert Z from weighted projective (Jacobian) to affine
    //i.e. (X, Y, Z) --> (X/Z^2, Y/Z^3)
    //also sets z to 1
    void to_affine(void)
    {
	element_invert(z, z);
	element_square(e0, z);
	element_mul(Z->x, Z->x, e0);
	element_mul(e0, e0, z);
	element_mul(Z->y, Z->y, e0);
	element_set1(z);
	element_set1(z2);
    }

    void proj_double(void)
    {
	const element_ptr x = Z->x;
	const element_ptr y = Z->y;
	//e0 = 3x^2 + (cc->a) z^4
	element_square(e0, x);
	element_mul_si(e0, e0, 3);
	element_square(e1, z2);
	element_mul(e1, e1, cca);
	element_add(e0, e0, e1);

	//z_out = 2 y z
	element_mul(z, y, z);
	element_mul_si(z, z, 2);
	element_square(z2, z);

	//e1 = 4 x y^2
	element_square(e2, y);
	element_mul(e1, x, e2);
	element_mul_si(e1, e1, 4);

	//x_out = e0^2 - 2 e1
	element_mul_si(e3, e1, 2);
	element_square(x, e0);
	element_sub(x, x, e3);

	//e2 = 8y^4
	element_square(e2, e2);
	element_mul_si(e2, e2, 8);

	//y_out = e0(e1 - x_out) - e2
	element_sub(e1, e1, x);
	element_mul(e0, e0, e1);
	element_sub(y, e0, e2);
    }

    void do_tangent(element_t e, element_t edenom, point_ptr A)
    {
	const element_ptr Zx = A->x;
	const element_ptr Zy = A->y;

	//a = -(3x^2 + cca z^4)
	//b = 2 y z^3
	//c = -(2 y^2 + x a)
	//a = z^2 a
	element_square(a, z2);
	element_mul(a, a, cca);
	element_square(b, Zx);
	element_mul_si(b, b, 3);
	element_add(a, a, b);
	element_neg(a, a);

	element_mul_si(e0, Zy, 2);
	element_mul(b, e0, z2);
	element_mul(b, b, z);

	element_mul(c, Zx, a);
	element_mul(a, a, z2);
	element_mul(e0, e0, Zy);
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

    void do_vertical(element_t e, element_t edenom, point_ptr A)
    {
	element_mul(e0, numx, z2);
	element_sub(e0, e0, A->x);
	element_mul(e, e, e0);

	element_mul(e0, denomx, z2);
	element_sub(e0, e0, A->x);
	element_mul(edenom, edenom, e0);
    }

    void do_line(element_ptr e, element_ptr edenom, point_ptr A, point_ptr B)
    {
	const element_ptr Ax = A->x;
	const element_ptr Ay = A->y;
	const element_ptr Bx = B->x;
	const element_ptr By = B->y;

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
    element_init(z, res->field);
    element_init(z2, res->field);
    element_set1(z);
    element_set1(z2);

    element_init(v, res->field);
    element_init(vd, res->field);
    element_init(v1, res->field);
    element_init(vd1, res->field);
    point_init(Z, P->curve);
    point_init(Z1, P->curve);

    point_set(Z, P);

    element_set1(v);
    element_set1(vd);
    element_set1(v1);
    element_set1(vd1);

    n = p->exp1;
    for (i=0; i<n; i++) {
	element_square(v, v);
	element_square(vd, vd);
	do_tangent(v, vd, Z);
	proj_double();
	do_vertical(vd, v, Z);
    }
    to_affine();
    if (p->sign1 < 0) {
	element_set(v1, vd);
	element_set(vd1, v);
	do_vertical(vd1, v1, Z);
	point_neg(Z1, Z);
    } else {
	element_set(v1, v);
	element_set(vd1, vd);
	point_set(Z1, Z);
    }
    n = p->exp2;
    for (; i<n; i++) {
	element_square(v, v);
	element_square(vd, vd);
	do_tangent(v, vd, Z);
	proj_double();
	do_vertical(vd, v, Z);
    }
    to_affine();
    element_mul(v, v, v1);
    element_mul(vd, vd, vd1);
    do_line(v, vd, Z, Z1);
    point_add(Z, Z, Z1);
    do_vertical(vd, v, Z);

    if (p->sign0 > 0) {
	do_vertical(v, vd, P);
    }

    element_invert(vd, vd);
    element_mul(res, v, vd);

    element_clear(v);
    element_clear(vd);
    element_clear(v1);
    element_clear(vd1);
    element_clear(z);
    element_clear(z2);
    point_clear(Z);
    point_clear(Z1);
    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(e0);
    element_clear(e1);
}

static void e_miller(element_t res, point_t P,
	element_ptr numx, element_ptr numy,
	element_ptr denomx, element_ptr denomy,
	e_pairing_data_ptr p)
{
    //collate divisions
    int n;
    element_t v, vd;
    element_t v1, vd1;
    point_t Z, Z1;
    element_t a, b, c;
    const common_curve_ptr cc = P->curve->data;
    element_t e0, e1;
    int i;

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
    element_init(v1, res->field);
    element_init(vd1, res->field);
    point_init(Z, P->curve);
    point_init(Z1, P->curve);

    point_set(Z, P);

    element_set1(v);
    element_set1(vd);
    element_set1(v1);
    element_set1(vd1);

    n = p->exp1;
    for (i=0; i<n; i++) {
	element_square(v, v);
	element_square(vd, vd);
	do_tangent(v, vd, Z);
	point_double(Z, Z);
	do_vertical(vd, v, Z);
    }
    if (p->sign1 < 0) {
	element_set(v1, vd);
	element_set(vd1, v);
	do_vertical(vd1, v1, Z);
	point_neg(Z1, Z);
    } else {
	element_set(v1, v);
	element_set(vd1, vd);
	point_set(Z1, Z);
    }
    n = p->exp2;
    for (; i<n; i++) {
	element_square(v, v);
	element_square(vd, vd);
	do_tangent(v, vd, Z);
	point_double(Z, Z);
	do_vertical(vd, v, Z);
    }
    element_mul(v, v, v1);
    element_mul(vd, vd, vd1);
    do_line(v, vd, Z, Z1);
    point_add(Z, Z, Z1);
    do_vertical(vd, v, Z);

    if (p->sign0 > 0) {
	do_vertical(v, vd, P);
    }

    element_invert(vd, vd);
    element_mul(res, v, vd);

    element_clear(v);
    element_clear(vd);
    element_clear(v1);
    element_clear(vd1);
    point_clear(Z);
    point_clear(Z1);
    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(e0);
    element_clear(e1);
}

static void e_pairing(element_ptr out, element_ptr in1, element_ptr in2,
	pairing_t pairing)
{
    e_pairing_data_ptr p = pairing->data;
    point_ptr Q = in2->data;
    point_t R, QR;
    point_init(R, p->Eq);
    point_init(QR, p->Eq);
    point_random(R);
    point_add(QR, Q, R);
    e_miller(out, in1->data, QR->x, QR->y, R->x, R->y, p);
    element_pow_mpz(out, out, p->tateexp);
    point_clear(R);
    point_clear(QR);
}

static void phi_identity(element_ptr out, element_ptr in, pairing_ptr pairing)
{
    (void) pairing;
    element_set(out, in);
}

void pairing_init_e_param(pairing_t pairing, e_param_t param)
{
    e_pairing_data_ptr p;
    element_t a, b;

    mpz_init(pairing->r);
    mpz_set(pairing->r, param->r);
    field_init_fp(pairing->Zr, pairing->r);
    pairing->map = e_pairing;

    p =	pairing->data = malloc(sizeof(e_pairing_data_t));
    p->exp2 = param->exp2;
    p->exp1 = param->exp1;
    p->sign1 = param->sign1;
    p->sign0 = param->sign0;
    field_init_fp(p->Fq, param->q);
    element_init(a, p->Fq);
    element_init(b, p->Fq);
    element_set_mpz(a, param->a);
    element_set_mpz(b, param->b);
    curve_init_cc_ab(p->Eq, a, b);

    mpz_init(p->tateexp);
    mpz_sub_ui(p->tateexp, p->Fq->order, 1);
    mpz_divexact(p->tateexp, p->tateexp, pairing->r);

    pairing->G1 = malloc(sizeof(field_t));
    pairing->G2 = pairing->G1;

    field_init_curve_group(pairing->G1, p->Eq, param->h);
    pairing->GT = p->Fq;
    pairing->phi = phi_identity;

    element_clear(a);
    element_clear(b);
}

//TODO: the following code is useless as the Tate pairing is degenerate on singular curves
static void sn_miller(element_t res, mpz_t q, point_t P,
	element_ptr Qx, element_ptr Qy)
{
    //collate divisions
    int m;
    element_t v, vd;
    point_t Z;
    element_t a, b, c;
    element_t e0, e1;

    void do_vertical(element_t e)
    {
	element_sub(e0, Qx, Z->x);
	element_mul(e, e, e0);
    }

    void do_tangent(element_t e)
    {
	//a = -slope_tangent(Z.x, Z.y);
	//b = 1;
	//c = -(Z.y + a * Z.x);
	//but we multiply by 2*Z.y to avoid division

	//a = -Zx * (Zx + Zx + Zx + 2)
	//b = 2 * Zy
	//c = -(2 Zy^2 + a Zx);
	element_ptr Zx = Z->x;
	element_ptr Zy = Z->y;

	element_mul_si(a, Zx, 3);
	element_set_si(e0, 2);
	element_add(a, a, e0);
	element_mul(a, a, Zx);
	element_neg(a, a);

	element_add(b, Zy, Zy);

	element_mul(e0, b, Zy);
	element_mul(c, a, Zx);
	element_add(c, c, e0);
	element_neg(c, c);

	element_mul(e0, a, Qx);
	element_mul(e1, b, Qy);
	element_add(e0, e0, e1);
	element_add(e0, e0, c);
	element_mul(e, e, e0);
    }

    void do_line(element_ptr e)
    {
	//a = -(B.y - A.y) / (B.x - A.x);
	//b = 1;
	//c = -(A.y + a * A.x);
	//but we'll multiply by B.x - A.x to avoid division

	element_ptr Ax = Z->x;
	element_ptr Ay = Z->y;
	element_ptr Bx = P->x;
	element_ptr By = P->y;

	element_sub(b, Bx, Ax);
	element_sub(a, Ay, By);
	element_mul(e0, b, Ay);
	element_mul(c, a, Ax);
	element_add(c, c, e0);
	element_neg(c, c);

	element_mul(e0, a, Qx);
	element_mul(e1, b, Qy);
	element_add(e0, e0, e1);
	element_add(e0, e0, c);
	element_mul(e, e, e0);
    }

    element_init(a, P->curve->field);
    element_init(b, P->curve->field);
    element_init(c, P->curve->field);
    element_init(e0, res->field);
    element_init(e1, res->field);

    element_init(v, res->field);
    element_init(vd, res->field);
    point_init(Z, P->curve);

    point_set(Z, P);

    element_set1(v);
    element_set1(vd);
    m = mpz_sizeinbase(q, 2) - 2;

    while(m >= 0) {
	element_mul(v, v, v);
	element_mul(vd, vd, vd);
	do_tangent(v);
	point_double(Z, Z);
	do_vertical(vd);
	if (mpz_tstbit(q, m)) {
	    do_line(v);
	    point_add(Z, Z, P);
	    do_vertical(vd);
	}
	m--;
    }

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

static void sn_pairing(element_ptr out, element_ptr in1, element_ptr in2,
	pairing_t pairing)
{
    e_pairing_data_ptr p = pairing->data;
    point_ptr Q = in2->data;
    element_t e0;
    point_t R, QR;
    point_init(R, p->Eq);
    point_init(QR, p->Eq);
    point_random(R);
    element_init(e0, out->field);
    point_add(QR, Q, R);
    sn_miller(out, pairing->r, in1->data, QR->x, QR->y);
    sn_miller(e0, pairing->r, in1->data, R->x, R->y);
    element_invert(e0, e0);
    element_mul(out, out, e0);
    //element_pow_mpz(out, out, p->tateexp);
    point_clear(R);
    point_clear(QR);
}

void pairing_init_singular_with_node(pairing_t pairing, mpz_t q)
{
    e_pairing_data_ptr p;
    mpz_t one;

    mpz_init(pairing->r);
    mpz_sub_ui(pairing->r, q, 1);
    field_init_fp(pairing->Zr, pairing->r);
    pairing->map = sn_pairing;

    p =	pairing->data = malloc(sizeof(e_pairing_data_t));
    field_init_fp(p->Fq, q);
    curve_init_singular_with_node(p->Eq, p->Fq);

    //mpz_init(p->tateexp);
    //mpz_sub_ui(p->tateexp, p->Fq->order, 1);
    //mpz_divexact(p->tateexp, p->tateexp, pairing->r);

    pairing->G1 = malloc(sizeof(field_t));
    pairing->G2 = pairing->G1;

    mpz_init(one);
    mpz_set_ui(one, 1);
    field_init_curve_group(pairing->G1, p->Eq, one);
    mpz_clear(one);
    pairing->GT = p->Fq;
    //pairing->phi = trace;
}
