#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>
#include "field.h"
#include "random.h"
#include "utils.h"
//Use Montgomery method for faster multiplication
//note inversion is slower
//an element of F_p is represented by xR (mod p)
//where R is the smallest power of the machine word greater than p
//TODO: why is this slower than naive implementation?

struct fp_field_data_s {
    size_t limbs;
    size_t bytes, bytes1;
    mp_limb_t *primelimbs;
    mp_limb_t *onelimbs;
    mp_limb_t pquote;
    //mp_limb_t *prelimbs;
};
typedef struct fp_field_data_s fp_field_data_t[1];
typedef struct fp_field_data_s *fp_field_data_ptr;

static void fp_init(element_ptr e)
{
    fp_field_data_ptr p = e->field->data;
    e->data = malloc(p->bytes);
}

static void fp_clear(element_ptr e)
{
    free(e->data);
}

/*
//Montgomery reduction, simple version
static inline void simple_mont_reduce(mp_limb_t *z, mp_limb_t *y, fp_field_data_ptr p)
{
    size_t t = p->limbs;
    mp_limb_t *tmp = malloc(p->bytes * 2);
    mp_limb_t *tmp1 = malloc(p->bytes * 2);
    mp_limb_t carry;
    mpn_mul_n(tmp, y, p->prelimbs, t);
    mpn_mul_n(tmp1, tmp, p->primelimbs, t);
    carry = mpn_add_n(tmp, tmp1, y, 2 * t);
    if (carry || mpn_cmp(&tmp[t], p->primelimbs, t) > 0) {
	mpn_sub_n(z, &tmp[t], p->primelimbs, t);
    } else {
	memcpy(z, &tmp[t], p->bytes);
    }
    free(tmp);
    free(tmp1);
}
*/

//Montgomery reduction. Modifies y
static inline void mont_reduce(mp_limb_t *z, mp_limb_t *y, fp_field_data_ptr p)
{
    size_t i;
    size_t t = p->limbs;
    //mp_limb_t *tmp = malloc(p->bytes1);
    mp_limb_t tmp[p->limbs];
    mp_limb_t carry;

    tmp[t] = mpn_mul_1(tmp, p->primelimbs, t, y[0] * p->pquote);
    //carry should increment at most once
    carry = mpn_add(&y[0], &y[0], 2 * t, tmp, t + 1);
    for (i=1; i<t; i++) {
	tmp[t] = mpn_mul_1(tmp, p->primelimbs, t, y[i] * p->pquote);
	//carry should increment at most once
	carry += mpn_add(&y[i], &y[i], 2 * t - i, tmp, t + 1);
    }
    if (carry || mpn_cmp(&y[t], p->primelimbs, t) > 0) {
	mpn_sub_n(z, &y[t], p->primelimbs, t);
    } else {
	memcpy(z, &y[t], p->bytes);
    }
    //free(tmp);
}

static inline void mont_reduce5(mp_limb_t *z, mp_limb_t *y, fp_field_data_ptr p)
{
    mp_limb_t tmp[6];
    mp_limb_t carry;

    tmp[5] = mpn_mul_1(tmp, p->primelimbs, 5, y[0] * p->pquote);
    //carry should increment at most once
    carry = mpn_add(&y[0], &y[0], 10, tmp, 6);

    tmp[5] = mpn_mul_1(tmp, p->primelimbs, 5, y[1] * p->pquote);
    carry += mpn_add(&y[1], &y[1], 2 * 5 - 1, tmp, 5 + 1);
    tmp[5] = mpn_mul_1(tmp, p->primelimbs, 5, y[2] * p->pquote);
    carry += mpn_add(&y[2], &y[2], 2 * 5 - 2, tmp, 5 + 1);
    tmp[5] = mpn_mul_1(tmp, p->primelimbs, 5, y[3] * p->pquote);
    carry += mpn_add(&y[3], &y[3], 2 * 5 - 3, tmp, 5 + 1);
    tmp[5] = mpn_mul_1(tmp, p->primelimbs, 5, y[4] * p->pquote);
    carry += mpn_add(&y[4], &y[4], 2 * 5 - 4, tmp, 5 + 1);

    if (carry || mpn_cmp(&y[5], p->primelimbs, 5) > 0) {
	mpn_sub_n(z, &y[5], p->primelimbs, 5);
    } else {
	memcpy(z, &y[5], 5 * sizeof(mp_limb_t));
    }
}

static void fp_to_mpz(mpz_ptr z, element_ptr e)
{
    fp_field_data_ptr p = e->field->data;
    mp_limb_t *tmp = malloc(p->bytes);
    mp_limb_t *y = malloc(p->bytes * 2);
    size_t t = p->limbs;
    memset(&y[t], 0, p->bytes);
    memcpy(y, e->data, p->bytes);
    mont_reduce(tmp, y, p);
    mpz_import(z, p->limbs, -1, sizeof(mp_limb_t), 0, 0, tmp);
    free(tmp);
    free(y);
}

static void fp_set_si(element_ptr e, signed long int op)
{
    if (op < 0) {
	printf("TODO: set_si negative!\n");
    }
    fp_field_data_ptr p = e->field->data;
    size_t n = p->limbs;
    mp_limb_t *d = e->data;
    mp_limb_t *tmp = malloc(sizeof(mp_limb_t) * (n + 1));
    mp_limb_t qp[2];
    memset(tmp, 0, sizeof(mp_limb_t) * n);
    tmp[n] = op;
    mpn_tdiv_qr(qp, d, 0, tmp, n + 1, p->primelimbs, n);
    free(tmp);
}

static void fp_set_mpz(element_ptr e, mpz_ptr z)
{
    if (mpz_sgn(z) < 0) {
	printf("TODO: set_mpz negative!\n");
    }
    size_t zlimbs = mpz_size(z);
    fp_field_data_ptr p = e->field->data;
    size_t n = p->limbs, check;
    mp_limb_t *d = e->data;
    mp_limb_t *tmp = malloc(sizeof(mp_limb_t) * (n + zlimbs));
    mp_limb_t *qp = malloc(sizeof(mp_limb_t) * (zlimbs + 1));
    memset(tmp, 0, sizeof(mp_limb_t) * n);
    mpz_export(&tmp[n], &check, -1, sizeof(mp_limb_t), 0, 0, z);
    //assert(check = zlimbs);
    mpn_tdiv_qr(qp, d, 0, tmp, n + zlimbs, p->primelimbs, n);
    free(qp);
    free(tmp);
}

static void fp_set0(element_ptr e)
{
    fp_field_data_ptr p = e->field->data;
    memset(e->data, 0, p->bytes);
}

static void fp_set1(element_ptr e)
{
    fp_field_data_ptr p = e->field->data;
    memcpy(e->data, p->onelimbs, p->bytes);
}

static size_t fp_out_str(FILE *stream, int base, element_ptr e)
{
    int res;
    mpz_t z;
    mpz_init(z);
    fp_to_mpz(z, e);
    res = mpz_out_str(stream, base, z);
    mpz_clear(z);
    return res;
}

static void fp_add(element_ptr r, element_ptr a, element_ptr b)
{
    fp_field_data_ptr p = r->field->data;
    size_t n = p->limbs;
    mp_limb_t carry;
    carry = mpn_add_n(r->data, a->data, b->data, n);

    if (carry || mpn_cmp(r->data, p->primelimbs, n) > 0) {
	mpn_sub_n(r->data, r->data, p->primelimbs, n);
    }
}

static void fp_sub(element_ptr r, element_ptr a, element_ptr b)
{
    fp_field_data_ptr p = r->field->data;
    size_t n = p->limbs;
    if (mpn_sub_n(r->data, a->data, b->data, n)) {
	mpn_add_n(r->data, r->data, p->primelimbs, n);
    }
}

static void fp_mul(element_ptr r, element_ptr a, element_ptr b)
{
    fp_field_data_ptr p = r->field->data;
    mp_limb_t *tmp = malloc(2 * p->bytes);
    mpn_mul_n(tmp, a->data, b->data, p->limbs);
    mont_reduce(r->data, tmp, p);
    free(tmp);
}

static void fp_mul5(element_ptr r, element_ptr a, element_ptr b)
{
    fp_field_data_ptr p = r->field->data;
    mp_limb_t *tmp = malloc(2 * p->bytes);
    mpn_mul_n(tmp, a->data, b->data, p->limbs);
    mont_reduce5(r->data, tmp, p);
    free(tmp);
}

/*
//really slow for some reason
static void fp_mul(element_ptr r, element_ptr a, element_ptr b)
{
    fp_field_data_ptr p = r->field->data;
    size_t i, t = p->limbs;
    mp_limb_t *tmp = malloc((4 * p->limbs + 4) * sizeof(mp_limb_t));
    mp_limb_t *tmp1 = &tmp[2 * p->limbs + 1];
    mp_limb_t *tmp2 = &tmp1[p->limbs + 2];
    mp_limb_t *z = r->data;
    mp_limb_t *x = a->data;
    mp_limb_t *y = b->data;
    mp_limb_t u;

    memset(tmp, 0, 2 * p->bytes + sizeof(mp_limb_t));

    u = (x[0] * y[0]) * p->pquote;
    tmp1[t] = mpn_mul_1(tmp1, y, t, x[0]);
    tmp2[t] = mpn_mul_1(tmp2, p->primelimbs, t, u);
    tmp1[t + 1] = mpn_add_n(tmp1, tmp1, tmp2, t + 1);
    mpn_add_n(&tmp[0], &tmp[0], tmp1, t + 2);

    for (i=1; i<t; i++) {
	u = (tmp[i] + x[i] * y[0]) * p->pquote;
	tmp1[t] = mpn_mul_1(tmp1, y, t, x[i]);
	tmp2[t] = mpn_mul_1(tmp2, p->primelimbs, t, u);
	tmp1[t + 1] = mpn_add_n(tmp1, tmp1, tmp2, t + 1);
	mpn_add_n(&tmp[i], &tmp[i], tmp1, t + 2);
    }
    if (tmp[2 * t] || mpn_cmp(&tmp[t], p->primelimbs, t) > 0) {
	mpn_sub_n(z, &tmp[t], p->primelimbs, t);
    } else {
	memcpy(z, &tmp[t], p->bytes);
    }
    free(tmp);
}
*/

static void fp_mul_mpz(element_ptr n, element_ptr a, mpz_ptr z)
{
    mpz_mul(n->data, a->data, z);
    mpz_mod(n->data, n->data, n->field->order);
}

static void fp_mul_si(element_ptr n, element_ptr a, signed long int z)
{
    mpz_mul_si(n->data, a->data, z);
    mpz_mod(n->data, n->data, n->field->order);
}

static void fp_set(element_ptr n, element_ptr a)
{
    mpz_set(n->data, a->data);
}

static void fp_neg(element_ptr n, element_ptr a)
{
    if (mpz_is0(a->data)) {
	mpz_set_ui(n->data, 0);
    } else {
	mpz_sub(n->data, n->field->order, a->data);
    }
}

static void fp_invert(element_ptr n, element_ptr a)
{
    UNUSED_VAR(n); UNUSED_VAR(a);
    /*
    fp_field_data_ptr p = a->field->data;
    mpz_invert(n->data, a->data, n->field->order);
    mpz_mul(n->data, n->data, p->R3modp);
    mont_reduce(n->data, n);
    */
}

static void fp_random(element_ptr n)
{
    //TODO: use random bytes directly?
    mpz_t z;
    mpz_init(z);
    pbc_mpz_random(z, n->field->order);
    fp_set_mpz(n, z);
    mpz_clear(z);
}

static void fp_from_hash(element_ptr n, int len, void *data)
    //TODO: something more sophisticated!
{
    mpz_t z;

    mpz_init(z);
    mpz_import(z, len, 1, 1, 0, 0, data);
    fp_set_mpz(n, z);
    mpz_clear(z);
}

static int fp_is1(element_ptr e)
{
    fp_field_data_ptr p = e->field->data;
    return !mpn_cmp(e->data, p->onelimbs, p->limbs);
}

static int fp_is0(element_ptr e)
{
    fp_field_data_ptr p = e->field->data;
    size_t i, t = p->limbs;
    mp_limb_t *d = e->data;
    for (i=0; i<t; i++) if (d[i]) return 0;
    return 1;
}

static int fp_cmp(element_ptr a, element_ptr b)
{
    fp_field_data_ptr p = a->field->data;
    size_t t = p->limbs;
    return mpn_cmp(a->data, b->data, t);
}

static int fp_is_sqr(element_ptr a)
{
    mpz_t z;
    int res;

    //0 is a square
    if (mpz_is0(a->data)) return 1;
    mpz_init(z);
    //mont_reduce(z, a);
    res = mpz_legendre(z, a->field->order) == 1;
    mpz_clear(z);
    return res;
}

static void fp_tonelli(element_ptr x, element_ptr a)
{
    int s;
    int i;
    mpz_t e;
    mpz_t t, t0;
    element_t ginv, e0;
    element_ptr nqr;

    mpz_init(t);
    mpz_init(e);
    mpz_init(t0);
    element_init(ginv, a->field);
    element_init(e0, a->field);
    nqr = field_get_nqr(a->field);

    element_invert(ginv, nqr); 

    //let q be the order of the field
    //q - 1 = 2^s t, t odd
    mpz_sub_ui(t, a->field->order, 1);
    s = mpz_scan1(t, 0);
    mpz_tdiv_q_2exp(t, t, s);
    mpz_set_ui(e, 0);
    for (i=2; i<=s; i++) {
	mpz_sub_ui(t0, a->field->order, 1);
	mpz_tdiv_q_2exp(t0, t0, i);
	element_pow_mpz(e0, ginv, e);
	element_mul(e0, e0, a);
	element_pow_mpz(e0, e0, t0);
	if (!element_is1(e0)) mpz_setbit(e, i-1);
    }
    element_pow_mpz(e0, ginv, e);
    element_mul(e0, e0, a);
    mpz_add_ui(t, t, 1);
    mpz_tdiv_q_2exp(t, t, 1);
    element_pow_mpz(e0, e0, t);
    mpz_tdiv_q_2exp(e, e, 1);
    element_pow_mpz(x, nqr, e);
    /* TODO: this would be a good place to use element_pow2 ... -hs */
    element_mul(x, x, e0);
    mpz_clear(t);
    mpz_clear(e);
    mpz_clear(t0);
    element_clear(ginv);
    element_clear(e0);
}

static void fp_field_clear(field_t f)
{
    fp_field_data_ptr p = f->data;
    free(p->primelimbs);
    free(p->onelimbs);
    free(p);
}

static int fp_to_bytes(unsigned char *data, element_t e)
{
    mpz_t z;
    int i, n;
    unsigned char *ptr;

    mpz_init(z);
    mpz_set(z, e->data);
    n = e->field->fixed_length_in_bytes;
    ptr = data;
    for (i = 0; i < n; i++) {
	*ptr = (unsigned char) mpz_get_ui(z);
	ptr++;
	mpz_tdiv_q_2exp(z, z, 8);
    }
    mpz_clear(z);
    return n;
}

static int fp_from_bytes(element_t e, unsigned char *data)
{
    unsigned char *ptr;
    int i, n;
    mpz_ptr z = e->data;
    mpz_t z1;

    mpz_init(z1);
    mpz_set_ui(z, 0);

    ptr = data;
    n = e->field->fixed_length_in_bytes;
    for (i=0; i<n; i++) {
	mpz_set_ui(z1, *ptr);
	mpz_mul_2exp(z1, z1, i * 8);
	ptr++;
	mpz_add(z, z, z1);
    }
    mpz_clear(z1);
    return n;
}

void field_init_mont_fp(field_ptr f, mpz_t prime)
{
    fp_field_data_ptr p;

    field_init(f);
    f->init = fp_init;
    f->clear = fp_clear;
    f->set_si = fp_set_si;
    f->set_mpz = fp_set_mpz;
    f->out_str = fp_out_str;
    f->add = fp_add;
    f->sub = fp_sub;
    f->set = fp_set;
    f->mul = fp_mul;
    f->mul_mpz = fp_mul_mpz;
    f->mul_si = fp_mul_si;
    f->neg = fp_neg;
    f->cmp = fp_cmp;
    f->invert = fp_invert;
    f->random = fp_random;
    f->from_hash = fp_from_hash;
    f->is1 = fp_is1;
    f->is0 = fp_is0;
    f->set0 = fp_set0;
    f->set1 = fp_set1;
    f->is_sqr = fp_is_sqr;
    f->sqrt = fp_tonelli;
    f->field_clear = fp_field_clear;
    f->to_bytes = fp_to_bytes;
    f->from_bytes = fp_from_bytes;
    f->to_mpz = fp_to_mpz;

    mpz_set(f->order, prime);
    f->fixed_length_in_bytes = (mpz_sizeinbase(prime, 2) + 7) / 8;
    p = f->data = malloc(sizeof(fp_field_data_t));
    p->limbs = mpz_size(prime);
    p->bytes = p->limbs * sizeof(mp_limb_t);
    p->bytes1 = (p->limbs + 1) * sizeof(mp_limb_t);
    p->primelimbs = malloc(p->bytes);
    p->onelimbs = malloc(p->bytes);
    mpz_export(p->primelimbs, &p->limbs, -1, sizeof(mp_limb_t), 0, 0, prime);

    if (p->limbs == 5) {
	f->mul = fp_mul5;
    }
    {
	size_t n = p->limbs;
	mp_limb_t *tmp = malloc(sizeof(mp_limb_t) * (n + 1));
	mp_limb_t qp[2];
	memset(tmp, 0, sizeof(mp_limb_t) * n);
	tmp[n] = 1;
	mpn_tdiv_qr(qp, p->onelimbs, 0, tmp, n + 1, p->primelimbs, n);
	free(tmp);
    }
    {
	mpz_t z, b;
	mpz_init(b);
	mpz_init(z);
	mpz_setbit(b, sizeof(mp_limb_t) * 8);
	mpz_invert(z, prime, b);
	mpz_sub(z, b, z);
	p->pquote = mpz_getlimbn(z, 0);

	/*
	mpz_set_ui(b, 0);
	mpz_setbit(b, sizeof(mp_limb_t) * 8 * p->limbs);
	mpz_invert(z, prime, b);
	mpz_sub(z, b, z);
	p->prelimbs = malloc(p->bytes);
	mpz_export(p->prelimbs, &p->limbs, -1, sizeof(mp_limb_t), 0, 0, z);
	*/

	mpz_clear(z);
	mpz_clear(b);
    }
}
