#include <alloca.h>
#include <string.h>
#include "field.h"
//Naive implementation of F_p
//using lowlevel GMP routines (mpn_* functions)
//Montgomery method seems to slow things down

struct fp_field_data_s {
    size_t limbs;
    size_t bytes;
    mp_limb_t *primelimbs;
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

static void fp_set_si(element_ptr e, signed long int op)
{
    if (op < 0) {
	printf("TODO: set_si negative!\n");
    }
    fp_field_data_ptr p = e->field->data;
    size_t t = p->limbs;
    mp_limb_t *d = e->data;
    d[0] = op;
    memset(&d[1], 0, sizeof(mp_limb_t) * (t - 1));
    //TODO have specialized functions for t = 1 case
    if (p->limbs == 1) d[0] = op % p->primelimbs[0];
}

static inline void from_mpz(element_ptr e, mpz_ptr z)
{
    fp_field_data_ptr p = e->field->data;
    size_t count;
    mpz_export(e->data, &count, -1, sizeof(mp_limb_t), 0, 0, z);
    memset(e->data + count * sizeof(mp_limb_t), 0, (p->limbs - count) * sizeof(mp_limb_t));
}

static void fp_to_mpz(mpz_ptr z, element_ptr a)
{
    fp_field_data_ptr p = a->field->data;
    mpz_import(z, p->limbs, -1, sizeof(mp_limb_t), 0, 0, a->data);
}

static void fp_set_mpz(element_ptr e, mpz_ptr z)
{
    mpz_t tmp;
    mpz_init(tmp);
    mpz_mod(tmp, z, e->field->order);
    from_mpz(e, tmp);
    mpz_clear(tmp);
}

static void fp_set0(element_ptr e)
{
    fp_field_data_ptr p = e->field->data;
    memset(e->data, 0, p->bytes);
}

static void fp_set1(element_ptr e)
{
    fp_field_data_ptr p = e->field->data;
    mp_limb_t *d = e->data;
    memset(d, 0, p->bytes);
    d[0] = 1;
}

static int fp_is1(element_ptr e)
{
    fp_field_data_ptr p = e->field->data;
    size_t i, t = p->limbs;
    mp_limb_t *d = e->data;
    for (i=1; i<t; i++) if (d[i]) return 0;
    return d[0] == 1;
}

static int fp_is0(element_ptr e)
{
    fp_field_data_ptr p = e->field->data;
    size_t i, t = p->limbs;
    mp_limb_t *d = e->data;
    for (i=0; i<t; i++) if (d[i]) return 0;
    return 1;
}

static size_t fp_out_str(FILE *stream, int base, element_ptr e)
{
    size_t result;
    mpz_t z;
    mpz_init(z);
    fp_to_mpz(z, e);
    result = mpz_out_str(stream, base, z);
    mpz_clear(z);
    return result;
}

static void fp_add(element_ptr r, element_ptr a, element_ptr b)
{
    fp_field_data_ptr p = r->field->data;
    size_t t = p->limbs;
    mp_limb_t carry;
    carry = mpn_add_n(r->data, a->data, b->data, t);

    if (carry || mpn_cmp(r->data, p->primelimbs, t) > 0) {
	mpn_sub_n(r->data, r->data, p->primelimbs, t);
    }
}

static void fp_sub(element_ptr r, element_ptr a, element_ptr b)
{
    fp_field_data_ptr p = r->field->data;
    size_t t = p->limbs;
    if (mpn_sub_n(r->data, a->data, b->data, t)) {
	mpn_add_n(r->data, r->data, p->primelimbs, t);
    }
}

static void fp_mul(element_ptr r, element_ptr a, element_ptr b)
{
    fp_field_data_ptr p = r->field->data;
    size_t t = p->limbs;
    //mp_limb_t tmp[3 * t + 1];
    //mp_limb_t *qp = &tmp[2 * t];
    mp_limb_t tmp[2 * t];
    mp_limb_t qp[t + 1];
    //mp_limb_t *tmp = alloca(p->bytes * 2);
    //mp_limb_t *qp = alloca(sizeof(mp_limb_t) * (t + 1));
    //mpn_mul_n(tmp, a->data, b->data, t);
    mpn_mul_n(tmp, a->data, b->data, t);

    mpn_tdiv_qr(qp, r->data, 0, tmp, 2 * t, p->primelimbs, t);
}

static void fp_mul_mpz(element_ptr e, element_ptr a, mpz_ptr op)
{
    mpz_t z;
    mpz_init(z);
    fp_to_mpz(z, a);
    mpz_mul(z, z, op);
    mpz_mod(z, z, e->field->order);
    from_mpz(e, z);

    mpz_clear(z);
}

static void fp_mul_si(element_ptr e, element_ptr a, signed long int op)
{
    if (op < 0) {
	printf("TODO: mul_si negative!\n");
    }
    fp_field_data_ptr p = e->field->data;
    size_t t = p->limbs;
    mp_limb_t tmp[t + 1];
    mp_limb_t qp[2];

    mpn_mul_1(tmp, a->data, t, op);
    mpn_tdiv_qr(qp, e->data, 0, tmp, t + 1, p->primelimbs, t);
}

static void fp_pow_mpz(element_ptr n, element_ptr a, mpz_ptr op)
{
    mpz_t z;
    mpz_init(z);
    fp_to_mpz(z, a);
    mpz_powm(z, z, op, n->field->order);
    from_mpz(n, op);
    mpz_clear(z);
}

static void fp_set(element_ptr n, element_ptr a)
{
    fp_field_data_ptr p = a->field->data;
    memcpy(n->data, a->data, p->bytes);
}

static void fp_neg(element_ptr n, element_ptr a)
{
    if (fp_is0(a->data)) {
	fp_set0(n->data);
    } else {
	fp_field_data_ptr p = a->field->data;
	mpn_sub_n(n->data, p->primelimbs, a->data, p->limbs);
    }
}

static void fp_invert(element_ptr e, element_ptr a)
{
    mpz_t z;
    mpz_init(z);
    fp_to_mpz(z, a);
    mpz_invert(z, z, e->field->order);
    from_mpz(e, z);
    mpz_clear(z);
}

static void fp_random(element_ptr n)
{
    mpz_t z;
    mpz_init(z);
    pbc_mpz_random(z, n->field->order);
    from_mpz(n, z);
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

static int fp_cmp(element_ptr a, element_ptr b)
{
    fp_field_data_ptr p = a->field->data;
    return mpn_cmp(a->data, b->data, p->limbs);
}

static int fp_is_sqr(element_ptr a)
{
    int res;
    mpz_t z;
    mpz_init(z);
    //0 is a square
    if (fp_is0(a)) return 1;
    fp_to_mpz(z, a);
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

static void fp_field_clear(field_t f)
{
    fp_field_data_ptr p = f->data;
    free(p->primelimbs);
    free(p);
}

void field_init_fast_fp(field_ptr f, mpz_t prime)
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
    f->pow_mpz = fp_pow_mpz;
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

    p = f->data = malloc(sizeof(fp_field_data_t));
    p->limbs = mpz_size(prime);
    p->bytes = p->limbs * sizeof(mp_limb_t);
    p->primelimbs = malloc(p->bytes);
    mpz_export(p->primelimbs, &p->limbs, -1, sizeof(mp_limb_t), 0, 0, prime);

    mpz_set(f->order, prime);
    f->fixed_length_in_bytes = (mpz_sizeinbase(prime, 2) + 7) / 8;
}
