#include <assert.h>
#include "a_param.h"
#include "param.h"
#include "tracker.h"

struct a_pairing_data_s {
    field_t Fq, Fq2;
    curve_t Eq, Eq2;
    mpz_t h;
    int exp2, exp1;
    int sign1;
};
typedef struct a_pairing_data_s a_pairing_data_t[1];
typedef struct a_pairing_data_s *a_pairing_data_ptr;

void a_param_init(a_param_t sp)
{
    mpz_init(sp->r);
    mpz_init(sp->q);
    mpz_init(sp->h);
}

void a_param_clear(a_param_t sp)
{
    mpz_clear(sp->r);
    mpz_clear(sp->q);
    mpz_clear(sp->h);
}

void a_param_gen(a_param_t sp, int rbits, int qbits)
{
    int found = 0;

    mpz_ptr q = sp->q;
    mpz_ptr r = sp->r;
    mpz_ptr h = sp->h;

    do {
	int i;
	mpz_set_ui(r, 0);

	if (rand() % 2) {
	    sp->exp2 = rbits - 1;
	    sp->sign1 = 1;
	} else {
	    sp->exp2 = rbits;
	    sp->sign1 = -1;
	}
	mpz_setbit(r, sp->exp2);

	//use q as a temp variable
	mpz_set_ui(q, 0);
	sp->exp1 = (rand() % (sp->exp2 - 1)) + 1;
	mpz_setbit(q, sp->exp1);
	if (sp->sign1 > 0) {
	    mpz_add(r, r, q);
	} else {
	    mpz_sub(r, r, q);
	}

	if (rand() % 2) {
	    sp->sign0 = 1;
	    mpz_add_ui(r, r, 1);
	} else {
	    sp->sign0 = -1;
	    mpz_sub_ui(r, r, 1);
	}
	if (!mpz_probab_prime_p(r, 10)) continue;
	for (i=0; i<10; i++) {
	    //use q as a temp variable
	    mpz_set_ui(q, 0);
	    mpz_setbit(q, qbits - rbits - 4 + 1);
	    pbc_mpz_random(h, q);
	    mpz_mul_ui(h, h, 12);
	    //finally q takes the value it should
	    mpz_mul(q, h, r);
	    mpz_sub_ui(q, q, 1);
	    if (mpz_probab_prime_p(q, 10)) {
		found = 1;
		break;
	    }
	}
    } while (!found);
}

void a_param_out_str(FILE *stream, a_param_ptr p)
{
    param_out_type(stream, "a");
    param_out_mpz(stream, "q", p->q);
    param_out_mpz(stream, "h", p->h);
    param_out_mpz(stream, "r", p->r);
    param_out_int(stream, "exp2", p->exp2);
    param_out_int(stream, "exp1", p->exp1);
    param_out_int(stream, "sign1", p->sign1);
    param_out_int(stream, "sign0", p->sign0);
}

void a_param_inp_generic (a_param_ptr p, fetch_ops_t *fops, void *ctx)
{
    assert (fops);
    assert (ctx);
    symtab_t tab;

    symtab_init(tab);
    param_read_generic (tab, fops, ctx);

    lookup_mpz(p->q, tab, "q");
    lookup_mpz(p->r, tab, "r");
    lookup_mpz(p->h, tab, "h");
    p->exp2 = lookup_int(tab, "exp2");
    p->exp1 = lookup_int(tab, "exp1");
    p->sign1 = lookup_int(tab, "sign1");
    p->sign0 = lookup_int(tab, "sign0");

    param_clear_tab(tab);
    symtab_clear(tab);
}

void a_param_inp_buf (a_param_ptr p, const char *buf, size_t len) 
{
    assert (buf);
    tracker_t t;
    tracker_init (&t, buf, len);
    a_param_inp_generic (p, &fops_buf, &t);
}

void a_param_inp_str (a_param_ptr p, FILE *stream)
{
    assert (stream);
    a_param_inp_generic (p, &fops_str, stream);
}

static void phi_identity(element_ptr out, element_ptr in, pairing_ptr pairing)
{
    (void) pairing;
    element_set(out, in);
}

static void a_pairing_proj(element_ptr out, element_ptr in1, element_ptr in2,
	pairing_t pairing)
//in1, in2 are from E(F_q), out from F_q^2
{
    a_pairing_data_ptr p = pairing->data;
    point_t V, V1;
    element_t z, z2;
    element_t f, f0, f1;
    element_t a, b, c;
    element_t e0;
    const element_ptr e1 = a, e2 = b, e3 = c;
    int i, n;
    element_ptr Qx = ((point_ptr) in2->data)->x;
    element_ptr Qy = ((point_ptr) in2->data)->y;

    //convert V from weighted projective (Jacobian) to affine
    //i.e. (X, Y, Z) --> (X/Z^2, Y/Z^3)
    //also sets z to 1
    void point_to_affine(void)
    {
	element_invert(z, z);
	element_square(e0, z);
	element_mul(V->x, V->x, e0);
	element_mul(e0, e0, z);
	element_mul(V->y, V->y, e0);
	element_set1(z);
	element_set1(z2);
    }

    void proj_double(void)
    {
	element_ptr x = V->x;
	element_ptr y = V->y;
	//e0 = 3x^2 + (cc->a) z^4
	element_square(e0, x);
	element_mul_si(e0, e0, 3);
	element_square(e1, z2);
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

    void do_tangent(void) {
	element_ptr Vx = V->x;
	element_ptr Vy = V->y;

	//a = -(3x^2 + cca z^4)
	//b = 2 y z^3
	//c = -(2 y^2 + x a)
	//a = z^2 a
	element_square(a, z2);
	element_square(b, Vx);
	element_mul_si(b, b, 3);
	element_add(a, a, b);
	element_neg(a, a);

	element_mul_si(e0, Vy, 2);
	element_mul(b, e0, z2);
	element_mul(b, b, z);

	element_mul(c, Vx, a);
	element_mul(a, a, z2);
	element_mul(e0, e0, Vy);
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

    void do_line(point_ptr A, point_ptr B) {
	//a = -(B.y - A.y) / (B.x - A.x);
	//b = 1;
	//c = -(A.y + a * A.x);
	//but we'll multiply by B.x - A.x to avoid division, so
	//a = -(By - Ay)
	//b = Bx - Ax
	//c = Ax By - Ay Bx;
	element_sub(a, A->y, B->y);
	element_sub(b, B->x, A->x);
	element_mul(c, A->x, B->y);
	element_mul(e0, A->y, B->x);
	element_sub(c, c, e0);

	//we'll map Q via (x,y) --> (-x, iy)
	//hence a Qx + c = -a Qx + c is real while
	//(b Qy) = b Qy i is purely imaginary.
	element_mul(a, a, Qx);
	element_sub(fi_re(f0), c, a);
	element_mul(fi_im(f0), b, Qy);
	element_mul(f, f, f0);
    }

    point_init(V, p->Eq);
    point_init(V1, p->Eq);
    point_set(V, in1->data);
    element_init(f, p->Fq2);
    element_init(f0, p->Fq2);
    element_init(f1, p->Fq2);
    element_set1(f);
    element_init(a, p->Fq);
    element_init(b, p->Fq);
    element_init(c, p->Fq);
    element_init(e0, p->Fq);
    element_init(z, p->Fq);
    element_init(z2, p->Fq);
    element_set1(z);
    element_set1(z2);
    n = p->exp1;
    for (i=0; i<n; i++) {
	//f = f^2 g_V,V(Q)
	//where g_V,V = tangent at V
	element_square(f, f);
	do_tangent();
	proj_double();
    }
    point_to_affine();
    if (p->sign1 < 0) {
	point_neg(V1, V);
	element_invert(f1, f);
    } else {
	point_set(V1, V);
	element_set(f1, f);
    }
    n = p->exp2;
    for (; i<n; i++) {
	element_square(f, f);
	do_tangent();
	proj_double();
    }

    element_mul(f, f, f1);
    point_to_affine();
    do_line(V, V1);

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
    element_clear(f1);
    element_clear(z);
    element_clear(z2);
    point_clear(V);
    point_clear(V1);
    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(e0);
}

static void a_pairing(element_ptr out, element_ptr in1, element_ptr in2,
	pairing_t pairing)
//in1, in2 are from E(F_q), out from F_q^2
{
    a_pairing_data_ptr p = pairing->data;
    point_t V, V1;
    element_t f, f0, f1;
    element_t a, b, c;
    element_t e0;
    int i, n;
    element_ptr Qx = ((point_ptr) in2->data)->x;
    element_ptr Qy = ((point_ptr) in2->data)->y;

    void do_tangent(void) {
	//a = -slope_tangent(V.x, V.y);
	//b = 1;
	//c = -(V.y + aV.x);
	//but we multiply by -2*V.y to avoid division so:
	//a = -(3 Vx^2 + cc->a)
	//b = 2 * Vy
	//c = -(2 Vy^2 + a Vx);
	element_ptr Vx = V->x;
	element_ptr Vy = V->y;
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

    void do_line(point_ptr A, point_ptr B) {
	//a = -(B.y - A.y) / (B.x - A.x);
	//b = 1;
	//c = -(A.y + a * A.x);
	//but we'll multiply by B.x - A.x to avoid division, so
	//a = -(By - Ay)
	//b = Bx - Ax
	//c = -(Ay b + a Ax);
	element_sub(a, A->y, B->y);
	element_sub(b, B->x, A->x);
	element_mul(c, A->x, B->y);
	element_mul(e0, A->y, B->x);
	element_sub(c, c, e0);

	//we'll map Q via (x,y) --> (-x, iy)
	//hence a Qx + c = -a Qx + c is real while
	//(b Qy) = b Qy i is purely imaginary.
	element_mul(a, a, Qx);
	element_sub(fi_re(f0), c, a);
	element_mul(fi_im(f0), b, Qy);
	element_mul(f, f, f0);
    }

    point_init(V, p->Eq);
    point_init(V1, p->Eq);
    point_set(V, in1->data);
    element_init(f, p->Fq2);
    element_init(f0, p->Fq2);
    element_init(f1, p->Fq2);
    element_set1(f);
    element_init(a, p->Fq);
    element_init(b, p->Fq);
    element_init(c, p->Fq);
    element_init(e0, p->Fq);
    n = p->exp1;
    element_t z; element_init(z, p->Fq);
    for (i=0; i<n; i++) {
	//f = f^2 g_V,V(Q)
	//where g_V,V = tangent at V
	element_square(f, f);
	do_tangent();
	point_double(V, V);
    }
    if (p->sign1 < 0) {
	point_neg(V1, V);
	element_invert(f1, f);
    } else {
	point_set(V1, V);
	element_set(f1, f);
    }
    n = p->exp2;
    for (; i<n; i++) {
	element_square(f, f);
	do_tangent();
	point_double(V, V);
    }

    element_mul(f, f, f1);
    do_line(V, V1);

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
    element_clear(f1);
    point_clear(V);
    point_clear(V1);
    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(e0);
}

static void pairing_clear_a_param(pairing_t pairing)
{
    a_pairing_data_ptr p = pairing->data;
    curve_clear(p->Eq);
    field_clear(p->Fq);
    field_clear(p->Fq2);
    mpz_clear(p->h);
    free(pairing->data);

    mpz_clear(pairing->r);
    field_clear(pairing->Zr);
    field_clear(pairing->G1);
    free(pairing->G1);
}

void pairing_init_a_param(pairing_t pairing, a_param_t param)
{
    element_t a, b;
    a_pairing_data_ptr p;

    p =	pairing->data = malloc(sizeof(a_pairing_data_t));
    p->exp2 = param->exp2;
    p->exp1 = param->exp1;
    p->sign1 = param->sign1;
    mpz_init(pairing->r);
    mpz_set(pairing->r, param->r);
    field_init_fp(pairing->Zr, pairing->r);
    pairing->map = a_pairing_proj;

    field_init_fp(p->Fq, param->q);
    element_init(a, p->Fq);
    element_init(b, p->Fq);
    element_set1(a);
    element_set0(b);
    curve_init_cc_ab(p->Eq, a, b);
    element_clear(a);
    element_clear(b);

    field_init_fi(p->Fq2, p->Fq);

    mpz_init(p->h);
    mpz_set(p->h, param->h);

    pairing->G1 = malloc(sizeof(field_t));
    field_init_curve_group(pairing->G1, p->Eq, param->h);
    pairing->G2 = pairing->G1;
    pairing->phi = phi_identity;
    pairing->GT = p->Fq2;
    pairing->clear_func = pairing_clear_a_param;
}
