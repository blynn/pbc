#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "field.h"
#include "darray.h"
#include "poly.h"
#include "curve.h"

size_t point_out_str(FILE *stream, int base, point_ptr p)
{
    size_t result, status;
    if (p->inf_flag) {
	if (EOF == fputc('O', stream)) return 0;
	return 1;
    }
    result = element_out_str(stream, base, p->x);
    if (!result) return 0;
    if (EOF == fputc(' ', stream)) return 0;
    status = element_out_str(stream, base, p->y);
    if (!status) return 0;
    return result + status + 1;
}

static void cc_from_x(point_ptr p, element_t x)
    //assumes there exists a point with given x coordinate
{
    element_t t;
    common_curve_ptr cc = p->curve->data;

    element_init(t, p->curve->field);
    p->inf_flag = 0;
    element_square(t, x);
    element_add(t, t, cc->a);
    element_mul(t, t, x);
    element_add(t, t, cc->b);
    element_sqrt(p->y, t);
    element_set(p->x, x);

    element_clear(t);
}

static void cc_random(point_ptr p)
{
    element_t t;
    common_curve_ptr cc = p->curve->data;

    element_init(t, p->curve->field);
    p->inf_flag = 0;
    do {
	element_random(p->x);
	element_square(t, p->x);
	element_add(t, t, cc->a);
	element_mul(t, t, p->x);
	element_add(t, t, cc->b);
    } while (!element_is_sqr(t));
    element_sqrt(p->y, t);

    element_clear(t);
}

static void cc_from_hash(point_ptr p, int len, void *data)
{
    //TODO: don't find a hash by the 255th try = freeze!
    void *datacopy;
    element_t t, t1;
    common_curve_ptr cc = p->curve->data;

    datacopy = malloc(len);
    memcpy(datacopy, data, len);

    element_init(t, p->curve->field);
    element_init(t1, p->curve->field);
    p->inf_flag = 0;
    for(;;) {
	element_from_hash(p->x, len, datacopy);
	element_square(t, p->x);
	element_add(t, t, cc->a);
	element_mul(t, t, p->x);
	element_add(t, t, cc->b);
	if (element_is_sqr(t)) break;
	((char *) datacopy)[0]++;
    }
    element_sqrt(p->y, t);

    element_clear(t);
    element_clear(t1);
    free(datacopy);
}

static void cc_neg(point_ptr r, point_ptr p)
{
    if (point_is_inf(p)) {
	point_set_inf(r);
	return;
    }
    point_set(r, p);
    element_neg(r->y, r->y);
}

static inline void cc_double_no_check(point_ptr r, point_ptr p)
{
    element_t lambda, e0, e1;
    common_curve_ptr cc = p->curve->data;

    element_init(lambda, p->curve->field);
    element_init(e0, p->curve->field);
    element_init(e1, p->curve->field);
    //same point: double them

    //lambda = (3x^2 + a) / 2y
    element_square(lambda, p->x);
    element_mul_si(lambda, lambda, 3);
    element_add(lambda, lambda, cc->a);
    element_add(e0, p->y, p->y);
    element_invert(e0, e0);
    element_mul(lambda, lambda, e0);
    //x1 = lambda^2 - 2x
    element_add(e1, p->x, p->x);
    element_square(e0, lambda);
    element_sub(e0, e0, e1);
    //y1 = (x - x1)lambda - y
    element_sub(e1, p->x, e0);
    element_mul(e1, e1, lambda);
    element_sub(e1, e1, p->y);

    element_set(r->x, e0);
    element_set(r->y, e1);
    r->inf_flag = 0;

    element_clear(lambda);
    element_clear(e0);
    element_clear(e1);
    return;
}

static inline void cc_double_no_check_ais0(point_ptr r, point_ptr p)
{
    element_t lambda, e0, e1;

    element_init(lambda, p->curve->field);
    element_init(e0, p->curve->field);
    element_init(e1, p->curve->field);
    //same point: double them

    //lambda = (3x^2 + a) / 2y
    element_square(lambda, p->x);
    element_mul_si(lambda, lambda, 3);
    element_add(e0, p->y, p->y);
    element_invert(e0, e0);
    element_mul(lambda, lambda, e0);
    //x1 = lambda^2 - 2x
    element_add(e1, p->x, p->x);
    element_square(e0, lambda);
    element_sub(e0, e0, e1);
    //y1 = (x - x1)lambda - y
    element_sub(e1, p->x, e0);
    element_mul(e1, e1, lambda);
    element_sub(e1, e1, p->y);

    element_set(r->x, e0);
    element_set(r->y, e1);
    r->inf_flag = 0;

    element_clear(lambda);
    element_clear(e0);
    element_clear(e1);
    return;
}

static void cc_double(point_ptr r, point_ptr p)
{
    if (point_is_inf(p)) {
	point_set_inf(r);
	return;
    }
    if (element_is0(p->y)) {
	point_set_inf(r);
	return;
    }
    //cc_double_no_check(r, p);
    r->curve->double_nocheck(r, p);
}

static void cc_add(point_ptr r, point_ptr p, point_ptr q)
{
    if (point_is_inf(p)) {
	point_set(r, q);
	return;
    }
    if (point_is_inf(q)) {
	point_set(r, p);
	return;
    }
    if (!element_cmp(p->x, q->x)) {
	if (!element_cmp(p->y, q->y)) {
	    if (element_is0(p->y)) {
		point_set_inf(r);
		return;
	    } else {
		//cc_double_no_check(r, p);
		r->curve->double_nocheck(r, p);
		return;
	    }
	}
	//points are inverses of each other
	point_set_inf(r);
	return;
    } else {
	element_t lambda, e0, e1;

	element_init(lambda, p->curve->field);
	element_init(e0, p->curve->field);
	element_init(e1, p->curve->field);

	//lambda = (y2-y1)/(x2-x1)
	element_sub(e0, q->x, p->x);
	element_invert(e0, e0);
	element_sub(lambda, q->y, p->y);
	element_mul(lambda, lambda, e0);
	//x3 = lambda^2 - x1 - x2
	element_square(e0, lambda);
	element_sub(e0, e0, p->x);
	element_sub(e0, e0, q->x);
	//y3 = (x1-x3)lambda - y1
	element_sub(e1, p->x, e0);
	element_mul(e1, e1, lambda);
	element_sub(e1, e1, p->y);

	element_set(r->x, e0);
	element_set(r->y, e1);
        r->inf_flag = 0;

	element_clear(lambda);
	element_clear(e0);
	element_clear(e1);
    }
}

static void cc_mul(point_ptr r, mpz_ptr n, point_ptr p)
{
    int s;

    point_t result;
    point_init(result, r->curve);
    point_set_inf(result);

    if (!mpz_is0(n)) for (s = mpz_sizeinbase(n, 2) - 1; s>=0; s--) {
	point_double(result, result);
	if (mpz_tstbit(n, s)) {
	    point_add(result, result, p);
	}
    }
    point_set(r, result);
    point_clear(result);
}

static void cc_clear(curve_ptr c)
{
    common_curve_ptr cc = c->data;
    element_clear(cc->a);
    element_clear(cc->b);
    free(c->data);
}

void cc_frobenius(point_ptr r, point_ptr p, mpz_ptr q)
{
    if (point_is_inf(p)) {
	point_set_inf(r);
	return;
    }
    element_pow_mpz(r->x, p->x, q);
    element_pow_mpz(r->y, p->y, q);
}

void curve_init_cc_ab(curve_ptr c, element_ptr a, element_ptr b)
{
    common_curve_ptr cc;
    c->field = a->field;
    c->random = cc_random;
    c->from_x = cc_from_x;
    c->from_hash = cc_from_hash;
    c->neg = cc_neg;
    if (element_is0(a)) {
	c->double_nocheck = cc_double_no_check_ais0;
    } else {
	c->double_nocheck = cc_double_no_check;
    }
    c->doublefn = cc_double;
    c->add = cc_add;
    c->mul = cc_mul;
    c->data = malloc(sizeof(common_curve_t));
    c->curve_clear = cc_clear;
    cc = c->data;
    element_init(cc->a, c->field);
    element_init(cc->b, c->field);
    element_set(cc->a, a);
    element_set(cc->b, b);
}

void curve_init_b(curve_ptr c, element_ptr b)
{
    element_t a;
    element_init(a, b->field);

    curve_init_cc_ab(c, a, b);

    element_clear(a);
}

void curve_init_cc_j(curve_ptr c, element_ptr j)
//assumes j != 0, 1728
{
    element_t a, b;
    element_init(a, j->field);
    element_init(b, j->field);

    element_set_si(a, 1728);
    element_sub(a, a, j);
    element_invert(a, a);
    element_mul(a, a, j);

    //b = 2 j / (1728 - j)
    element_add(b, a, a);
    //a = 3 j / (1728 - j)
    element_add(a, a, b);
    curve_init_cc_ab(c, a, b);

    element_clear(a);
    element_clear(b);
}

void twist_curve(curve_ptr c)
    //WARNING: existing points will no longer lie on c
    //as this modifies c in place
{
    common_curve_ptr cc = c->data;
    element_ptr nqr = field_get_nqr(c->field);
    element_mul(cc->a, cc->a, nqr);
    element_mul(cc->a, cc->a, nqr);
    element_mul(cc->b, cc->b, nqr);
    element_mul(cc->b, cc->b, nqr);
    element_mul(cc->b, cc->b, nqr);
}

void cc_init_map_curve(curve_ptr cnew, curve_ptr c,
	field_ptr dstfield, fieldmap map)
{
    common_curve_ptr ccnew, cc;
    cnew->field = dstfield;
    cnew->random = c->random;
    cnew->neg = c->neg;
    cnew->add = c->add;
    cnew->double_nocheck = c->double_nocheck;
    cnew->doublefn = c->doublefn;
    cnew->mul = c->mul;
    cnew->from_x = c->from_x;
    cnew->curve_clear = c->curve_clear;

    cnew->data = malloc(sizeof(common_curve_t));
    ccnew = cnew->data;
    cc = c->data;
    element_init(ccnew->a, cnew->field);
    element_init(ccnew->b, cnew->field);
    map(ccnew->a, cc->a);
    map(ccnew->b, cc->b);
}

void compute_trace_n(mpz_t res, mpz_t q, mpz_t trace, int n)
    //compute trace of Frobenius at q^n given trace at q
    //see p.105 of Blake, Seroussi and Smart
{
    int i;
    mpz_t c0, c1, c2;
    mpz_t t0;

    mpz_init(c0);
    mpz_init(c1);
    mpz_init(c2);
    mpz_init(t0);
    mpz_set_ui(c2, 2);
    mpz_set(c1, trace);
    for (i=2; i<=n; i++) {
	mpz_mul(c0, trace, c1);
	mpz_mul(t0, q, c2);
	mpz_sub(c0, c0, t0);
	mpz_set(c2, c1);
	mpz_set(c1, c0);
    }
    mpz_set(res, c1);
    mpz_clear(t0);
    mpz_clear(c2);
    mpz_clear(c1);
    mpz_clear(c0);
}

void point_map(point_t R, fieldmap map, point_t P)
{
    if (point_is_inf(P)) {
	point_set_inf(R);
	return;
    }
    map(R->x, P->x);
    map(R->y, P->y);
    R->inf_flag = 0;
}

static void curve_group_init(element_ptr e)
{
    e->data = malloc(sizeof(point_t));
    point_init(e->data, ((curve_group_data_ptr) e->field->data)->curve);
}

static void curve_group_clear(element_ptr e)
{
    point_clear(e->data);
    free(e->data);
}

static void curve_group_invert(element_ptr x, element_ptr a)
{
    point_neg(x->data, a->data);
}

static void curve_group_mul(element_ptr x, element_ptr a, element_ptr b)
{
    point_add(x->data, a->data, b->data);
}

static void curve_group_div(element_ptr x, element_ptr a, element_ptr b)
{
    if (a == b) {
	point_set_inf(x->data);
    } else {
	point_neg(b->data, b->data);
	point_add(x->data, a->data, b->data);
	point_neg(b->data, b->data);
    }
}

static void curve_group_set1(element_ptr x)
{
    point_set_inf(x->data);
}

static void curve_group_set(element_ptr x, element_ptr y)
{
    point_set(x->data, y->data);
}

static void curve_group_random(element_ptr x)
{
    curve_group_data_ptr p = x->field->data;
    point_random(x->data);
    point_mul(x->data, p->cofac, x->data);
}

static void curve_group_from_hash(element_ptr x, int len, void *data)
{
    curve_group_data_ptr p = x->field->data;
    point_from_hash(x->data, len, data);
    point_mul(x->data, p->cofac, x->data);
}

static size_t curve_group_out_str(FILE *stream, int base, element_ptr x)
{
    return point_out_str(stream, base, x->data);
}

static int curve_group_length_in_bytes(element_ptr x)
{
    point_ptr P = x->data;
    return element_length_in_bytes(P->x) + element_length_in_bytes(P->y);
}

static int curve_group_to_bytes(unsigned char *data, element_t e)
{
    point_ptr P = e->data;
    int len;
    len = element_to_bytes(data, P->x);
    len += element_to_bytes(data + len, P->y);
    return len;
}

static int curve_group_from_bytes(element_t e, unsigned char *data)
{
    point_ptr P = e->data;
    int len;

    //TODO: use point_set() instead of messing with internal data structure
    P->inf_flag = 0;
    len = element_from_bytes(P->x, data);
    len += element_from_bytes(P->y, data + len);
    return len;
}

static void curve_group_print_info(FILE *out, field_t f)
{
    int len;
    fprintf(out, "Group of points on elliptic curve");
    if ((len = f->fixed_length_in_bytes)) {
	fprintf(out, ", bits per coord = %d", len * 8 / 2);
    }
    fprintf(out, "\n");
}

void field_clear_curve_group(field_t f)
{
    curve_group_data_ptr p;
    p = f->data;
    mpz_clear(p->cofac);
    free(p);
}

void field_init_curve_group(field_t f, curve_t c, mpz_t cofac)
{
    curve_group_data_ptr p;
    field_init(f);
    p = f->data = malloc(sizeof(curve_group_data_t));
    p->curve = c;
    mpz_init(p->cofac);
    mpz_set(p->cofac, cofac);
    f->init = curve_group_init;
    f->clear = curve_group_clear;
    f->neg = f->invert = curve_group_invert;
    f->add = f->mul = curve_group_mul;
    f->sub = curve_group_div;
    f->set0 = f->set1 = curve_group_set1;
    f->set = curve_group_set;
    f->random = curve_group_random;
    f->from_hash = curve_group_from_hash;
    f->out_str = curve_group_out_str;
    f->field_clear = field_clear_curve_group;
    if (c->field->fixed_length_in_bytes < 0) {
	f->length_in_bytes = curve_group_length_in_bytes;
    } else {
	f->fixed_length_in_bytes = 2 * c->field->fixed_length_in_bytes;
    }
    f->to_bytes = curve_group_to_bytes;
    f->from_bytes = curve_group_from_bytes;
    f->print_info = curve_group_print_info;
}

//singular with node: y^2 = x^3 + x^2
static void sn_random(point_ptr p)
{
    element_t t;

    element_init(t, p->curve->field);
    p->inf_flag = 0;
    do {
	element_random(p->x);
	if (element_is0(p->x)) continue;
	element_square(t, p->x);
	element_add(t, t, p->x);
	element_mul(t, t, p->x);
    } while (!element_is_sqr(t));
    element_sqrt(p->y, t);

    element_clear(t);
}

static inline void sn_double_no_check(point_ptr r, point_ptr p)
{
    element_t lambda, e0, e1;

    element_init(lambda, p->curve->field);
    element_init(e0, p->curve->field);
    element_init(e1, p->curve->field);
    //same point: double them

    //lambda = (3x^2 + 2x) / 2y
    element_mul_si(lambda, p->x, 3);
    element_set_si(e0, 2);
    element_add(lambda, lambda, e0);
    element_mul(lambda, lambda, p->x);
    element_add(e0, p->y, p->y);
    element_invert(e0, e0);
    element_mul(lambda, lambda, e0);
    //x1 = lambda^2 - 2x - 1
    element_add(e1, p->x, p->x);
    element_square(e0, lambda);
    element_sub(e0, e0, e1);
    element_set_si(e1, 1);
    element_sub(e0, e0, e1);
    //y1 = (x - x1)lambda - y
    element_sub(e1, p->x, e0);
    element_mul(e1, e1, lambda);
    element_sub(e1, e1, p->y);

    element_set(r->x, e0);
    element_set(r->y, e1);
    r->inf_flag = 0;

    element_clear(lambda);
    element_clear(e0);
    element_clear(e1);
    return;
}

static void sn_double(point_ptr r, point_ptr p)
{
    if (point_is_inf(p)) {
	point_set_inf(r);
	return;
    }
    if (element_is0(p->y)) {
	point_set_inf(r);
	return;
    }
    sn_double_no_check(r, p);
}

static void sn_add(point_ptr r, point_ptr p, point_ptr q)
{
    if (point_is_inf(p)) {
	point_set(r, q);
	return;
    }
    if (point_is_inf(q)) {
	point_set(r, p);
	return;
    }
    if (!element_cmp(p->x, q->x)) {
	if (!element_cmp(p->y, q->y)) {
	    if (element_is0(p->y)) {
		point_set_inf(r);
		return;
	    } else {
		sn_double_no_check(r, p);
		return;
	    }
	}
	//points are inverses of each other
	point_set_inf(r);
	return;
    } else {
	element_t lambda, e0, e1;

	element_init(lambda, p->curve->field);
	element_init(e0, p->curve->field);
	element_init(e1, p->curve->field);

	//lambda = (y2-y1)/(x2-x1)
	element_sub(e0, q->x, p->x);
	element_invert(e0, e0);
	element_sub(lambda, q->y, p->y);
	element_mul(lambda, lambda, e0);
	//x3 = lambda^2 - x1 - x2 - 1
	element_square(e0, lambda);
	element_sub(e0, e0, p->x);
	element_sub(e0, e0, q->x);
	element_set1(e1);
	element_sub(e0, e0, e1);
	//y3 = (x1-x3)lambda - y1
	element_sub(e1, p->x, e0);
	element_mul(e1, e1, lambda);
	element_sub(e1, e1, p->y);

	element_set(r->x, e0);
	element_set(r->y, e1);
        r->inf_flag = 0;

	element_clear(lambda);
	element_clear(e0);
	element_clear(e1);
    }
}

static void sn_clear(curve_ptr c)
{
    (void) c;
}

//nonsingular points on sn curves map to finite field elements via
//(x, y) --> (y + x)/(y - x)
//and the reverse map is
//a --> (4a/(a-1)^2, 4a(a+1)/(a-1)^3)
void sn_point_to_field(element_t out, point_ptr P)
{
    element_t e0, e1;
    if (P->inf_flag) {
	element_set1(out);
	return;
    }
    element_init(e0, out->field);
    element_init(e1, out->field);
    element_add(e0, P->y, P->x);
    element_sub(e1, P->y, P->x);
    element_invert(e1, e1);
    element_mul(out, e0, e1);
    element_clear(e0);
    element_clear(e1);
}

void sn_field_to_point(point_ptr P, element_t in)
{
    element_t e0, e1, e2;

    if (element_is1(in)) {
	point_set_inf(P);
	return;
    }
    element_init(e0, in->field);
    element_init(e1, in->field);
    element_init(e2, in->field);

    element_set1(e1);
    element_sub(e0, in, e1);
    element_invert(e0, e0);

    element_mul_si(e2, in, 4);

    element_add(P->y, in, e1);

    element_mul(e1, e0, e0);
    element_mul(P->x, e1, e2);
    element_mul(P->y, P->y, e2);
    element_mul(P->y, P->y, e0);
    element_mul(P->y, P->y, e1);
    P->inf_flag = 0;

    element_clear(e0);
    element_clear(e1);
    element_clear(e2);
}

void curve_init_singular_with_node(curve_ptr c, field_t field)
{
    c->field = field;
    c->random = sn_random;
    //c->from_x = cc_from_x;
    //c->from_hash = cc_from_hash;
    c->neg = cc_neg;
    c->doublefn = sn_double;
    c->add = sn_add;
    c->mul = cc_mul;
    c->curve_clear = sn_clear;
}

int element_to_bytes_compressed(unsigned char *data, element_ptr e)
    //e must be a point on an elliptic curve
{
    point_ptr P = e->data;
    int len;
    len = element_to_bytes(data, P->x);
    if (element_sign(P->y) > 0) {
	data[len] = 1;
    } else {
	data[len] = 0;
    }
    len++;
    return len;
}

int element_from_bytes_compressed(element_ptr e, unsigned char *data)
    //e must be a point on an elliptic curve
{
    point_ptr P = e->data;
    int len;
    len = element_from_bytes(P->x, data);
    point_from_x(P, P->x);

    if (data[len]) {
	if (element_sign(P->y) < 0) element_neg(P->y, P->y);
    } else if (element_sign(P->y) > 0) {
	element_neg(P->y, P->y);
    }
    len++;
    return len;
}

int element_length_in_bytes_compressed(element_ptr e)
{
    point_ptr P = e->data;
    return element_length_in_bytes(P->x) + 1;
}

int element_to_bytes_x_only(unsigned char *data, element_ptr e)
    //e must be a point on an elliptic curve
{
    point_ptr P = e->data;
    int len;
    len = element_to_bytes(data, P->x);
    return len;
}

int element_from_bytes_x_only(element_ptr e, unsigned char *data)
    //e must be a point on an elliptic curve
{
    point_ptr P = e->data;
    int len;
    len = element_from_bytes(P->x, data);
    point_from_x(P, P->x);
    return len;
}

int element_length_in_bytes_x_only(element_ptr e)
{
    point_ptr P = e->data;
    return element_length_in_bytes(P->x);
}
