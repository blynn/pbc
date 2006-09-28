#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <alloca.h>
#include <string.h>
#include <gmp.h>
#include "field.h"
#include "random.h"
#include "fp.h"
// F_p using Montgomery representation
// Let b = 256^sizeof(mp_limb_t)
// Let R = b^t be the smallest power of b greater than the modulus p
// Then x is stored as xR (mod p)
// Additive operations are identical
// Multipicative ones use Montgomery reduction
// only works with odd modulus

struct fp_field_data_s {
    size_t limbs;
    size_t bytes;
    mp_limb_t *primelimbs;
    mp_limb_t negpinv;  //-p^-1 mod b
    mp_limb_t *R; //R mod p
};
typedef struct fp_field_data_s fp_field_data_t[1];
typedef struct fp_field_data_s *fp_field_data_ptr;

struct data_s {
    int flag;
    mp_limb_t *d;
};
typedef struct data_s *dataptr;

static void mont_reduce(mp_limb_t *x, mp_limb_t *y, fp_field_data_ptr p)
{
    //Algorithm II.4 from Blake, Seroussi and Smart
    size_t t = p->limbs;
    size_t i;
    mp_limb_t flag = 0;
    for (i=0; i<t; i++) {
	mp_limb_t u = y[i] * p->negpinv;
	mp_limb_t carry = mpn_addmul_1(&y[i], p->primelimbs, t, u);
	//mpn_add_1(&y[i+t], &y[i+t], t - i + 1, carry);
	flag += mpn_add_1(&y[i+t], &y[i+t], t - i, carry);
    }
    //TODO: do I need to check flag?
    if (flag || mpn_cmp(&y[t], p->primelimbs, t) >= 0) {
	mpn_sub_n(x, &y[t], p->primelimbs, t);
    } else {
	//TODO: GMP set
	memcpy(x, &y[t], t * sizeof(mp_limb_t));
    }
}

//store x as xR mod prime
//assumes z != 0
static void from_mpz(element_t e, mpz_t z)
{
    fp_field_data_ptr p = e->field->data;
    dataptr ed = e->data;
    size_t count;
    mpz_t z1;
    mpz_init(z1);
    mpz_mul_2exp(z1, z, p->bytes * 8);
    mpz_mod(z1, z1, e->field->order);

    mpz_export(ed->d, &count, -1, sizeof(mp_limb_t), 0, 0, z1);
    memset((void *) (((unsigned int) ed->d) + count * sizeof(mp_limb_t)),
	0, (p->limbs - count) * sizeof(mp_limb_t));
    mpz_clear(z1);
}

static void fp_init(element_ptr e)
{
    fp_field_data_ptr p = e->field->data;
    dataptr dp = e->data = malloc(sizeof(struct data_s));
    dp->flag = 0;
    dp->d = malloc(p->bytes);
}

static void fp_clear(element_ptr e)
{
    dataptr dp = e->data;
    free(dp->d);
    free(e->data);
}

static void fp_set_mpz(element_ptr e, mpz_ptr z)
{
    dataptr dp = e->data;
    if (!mpz_sgn(z)) {
	dp->flag = 0;
    } else {
	mpz_t tmp;
	mpz_init(tmp);
	mpz_mod(tmp, z, e->field->order);
	from_mpz(e, tmp);
	mpz_clear(tmp);
	dp->flag = 2;
    }
}

static void fp_set_si(element_ptr e, signed long int op)
{
    dataptr dp = e->data;
    if (!op) {
	dp->flag = 0;
    } else {
	const fp_field_data_ptr p = e->field->data;
	const size_t t = p->limbs;
	if (op < 0) {
	    mpn_sub_1(dp->d, p->primelimbs, t, -op);
	} else {
	    dp->d[0] = op;
	    memset(&dp->d[1], 0, sizeof(mp_limb_t) * (t - 1));
	}
	dp->flag = 2;
    }
}

//x is stored as xR, hence apply reduction before output
static void fp_to_mpz(mpz_ptr z, element_ptr e)
{
    dataptr dp = e->data;
    if (!dp->flag) {
	mpz_set_ui(z, 0);
    } else {
	fp_field_data_ptr p = e->field->data;
	mp_limb_t tmp[2 * p->limbs];

	memcpy(tmp, dp->d, p->limbs * sizeof(mp_limb_t));
	memset(&tmp[p->limbs], 0, p->limbs * sizeof(mp_limb_t));
	_mpz_realloc(z, p->limbs);
	mont_reduce(z->_mp_d, tmp, p);
	z->_mp_size = p->limbs;
	while(!z->_mp_d[z->_mp_size - 1]) {
	    z->_mp_size--;
	}
	//mpz_import(z, p->limbs, -1, sizeof(mp_limb_t), 0, 0, dp->d);
    }
}

static void fp_set0(element_ptr e)
{
    dataptr dp = e->data;
    dp->flag = 0;
}

static void fp_set1(element_ptr e)
{
    fp_field_data_ptr p = e->field->data;
    dataptr dp = e->data;
    dp->flag = 2;
    memcpy(dp->d, p->R, p->bytes);
}

static int fp_is1(element_ptr e)
{
    dataptr dp = e->data;
    if (!dp->flag) return 0;
    else {
	fp_field_data_ptr p = e->field->data;
	size_t i, t = p->limbs;
	if (dp->d[0] != 1) return 0;
	for (i=1; i<t; i++) if (dp->d[i]) return 0;
	return 1;
    }
}

static int fp_is0(element_ptr e)
{
    dataptr dp = e->data;
    return !dp->flag;
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

static void fp_set(element_ptr c, element_ptr a)
{
    dataptr ad = a->data;
    dataptr cd = c->data;
    if (a == c) return;
    if (!ad->flag) {
	cd->flag = 0;
    } else {
	fp_field_data_ptr p = a->field->data;

	//Assembly is faster here, but I don't want to stoop to that level.
	//Instead of calling slower memcpy, wrap stuff so that GMP assembly
	//gets called.
	/*
	memcpy(cd->d, ad->d, p->bytes);
	*/
	mpz_t z1, z2;
	z1->_mp_d = cd->d;
	z2->_mp_d = ad->d;
	z1->_mp_size = z1->_mp_alloc = z2->_mp_size = z2->_mp_alloc = p->limbs;
	mpz_set(z1, z2);

	cd->flag = 2;
    }
}

static void fp_add(element_ptr c, element_ptr a, element_ptr b)
{
    dataptr ad = a->data, bd = b->data;

    if (!ad->flag) {
	fp_set(c, b);
    } else if (!bd->flag) {
	fp_set(c, a);
    } else {
	dataptr cd = c->data;
	fp_field_data_ptr p = a->field->data;
	const size_t t = p->limbs;
	mp_limb_t carry;
	carry = mpn_add_n(cd->d, ad->d, bd->d, t);

	if (carry) {
	    //assumes result of following sub is not zero,
	    //i.e. modulus cannot be 2^(n * bits_per_limb)
	    mpn_sub_n(cd->d, cd->d, p->primelimbs, t);
	    cd->flag = 2;
	} else {
	    int i = mpn_cmp(cd->d, p->primelimbs, t);
	    if (!i) {
		cd->flag = 0;
	    } else {
		cd->flag = 2;
		if (i > 0) {
		    mpn_sub_n(cd->d, cd->d, p->primelimbs, t);
		}
	    }
	}
    }
}

static void fp_double(element_ptr c, element_ptr a)
{
    dataptr ad = a->data, cd = c->data;
    if (!ad->flag) {
	cd->flag = 0;
    } else {
	fp_field_data_ptr p = c->field->data;
	const size_t t = p->limbs;
	if (mpn_lshift(cd->d, ad->d, t, 1)) {
	    cd->flag = 2;
	    //again, assumes result is not zero:
	    mpn_sub_n(cd->d, cd->d, p->primelimbs, t);
	} else  {
	    int i = mpn_cmp(cd->d, p->primelimbs, t);
	    if (!i) {
		cd->flag =0;
	    } else {
		cd->flag = 2;
		if (i > 0) {
		    mpn_sub_n(cd->d, cd->d, p->primelimbs, t);
		}
	    }
	}
    }
}

static void fp_neg(element_ptr c, element_ptr a)
{
    dataptr ad = a->data, cd = c->data;
    if (!ad->flag) cd->flag = 0;
    else {
	fp_field_data_ptr p = a->field->data;
	mpn_sub_n(cd->d, p->primelimbs, ad->d, p->limbs);
	cd->flag = 2;
    }
}

static void fp_sub(element_ptr c, element_ptr a, element_ptr b)
{
    dataptr ad = a->data, bd = b->data;

    if (!ad->flag) {
	fp_neg(c, b);
    } else if (!bd->flag) {
	fp_set(c, a);
    } else {
	fp_field_data_ptr p = c->field->data;
	size_t t = p->limbs;
	dataptr cd = c->data;
	int i = mpn_cmp(ad->d, bd->d, t);

	if (i == 0) {
	    cd->flag = 0;
	} else {
	    cd->flag = 2;
	    mpn_sub_n(cd->d, ad->d, bd->d, t);
	    if (i < 0) {
		mpn_add_n(cd->d, cd->d, p->primelimbs, t);
	    }
	}
    }
}

static void fp_mul(element_ptr c, element_ptr a, element_ptr b)
{
    dataptr ad = a->data, bd = b->data;
    dataptr cd = c->data;

    if (!ad->flag || !bd->flag) {
	cd->flag = 0;
    } else {
	if (0) {
	//TODO: z[t] could overflow?
	fp_field_data_ptr p = c->field->data;
	size_t i, t = p->limbs;
	mp_limb_t z[t+1];
	mp_limb_t u = (ad->d[0] * bd->d[0]) * p->negpinv;
	z[t] = mpn_mul_1(z, bd->d, t, ad->d[0]);
//printf("d0: %lu\n", z[t]);
	z[t] += mpn_addmul_1(z, p->primelimbs, t, u);
//printf("d0: %lu\n", z[t]);
	memmove(z, &z[1], t * sizeof(mp_limb_t));
	for (i=1; i<t; i++) {
	    u = (z[0] + ad->d[i] * bd->d[0]) * p->negpinv;
	    z[t] = mpn_addmul_1(z, bd->d, t, ad->d[i]);
//printf("d0: %lu\n", z[t]);
	    z[t] += mpn_addmul_1(z, p->primelimbs, t, u);
//printf("d0: %lu\n", z[t]);
	    memmove(z, &z[1], t * sizeof(mp_limb_t));
	}
	if (mpn_cmp(z, p->primelimbs, t) >= 0) {
	    mpn_sub_n(cd->d, z, p->primelimbs, t);
	} else {
	    //TODO: GMP set
	    memcpy(cd->d, z, t * sizeof(mp_limb_t));
	}
	}

	//TODO: z[t] could overflow?
	fp_field_data_ptr p = c->field->data;
	size_t i, t = p->limbs;
	mp_limb_t z[2 * t+1];
	mp_limb_t u = (ad->d[0] * bd->d[0]) * p->negpinv;
	z[t] = mpn_mul_1(z, bd->d, t, ad->d[0]);
//printf("d0: %lu\n", z[t]);
	z[t] += mpn_addmul_1(z, p->primelimbs, t, u);
//printf("d0: %lu\n", z[t]);
	for (i=1; i<t; i++) {
	    u = (z[i] + ad->d[i] * bd->d[0]) * p->negpinv;
	    z[t + i] = mpn_addmul_1(z + i, bd->d, t, ad->d[i]);
//printf("d0: %lu\n", z[t]);
	    z[t + i] += mpn_addmul_1(z + i, p->primelimbs, t, u);
//printf("d0: %lu\n", z[t]);
	}
	if (mpn_cmp(z + t, p->primelimbs, t) >= 0) {
	    mpn_sub_n(cd->d, z + t, p->primelimbs, t);
	} else {
	    //TODO: GMP set
	    memcpy(cd->d, z + t, t * sizeof(mp_limb_t));
	}

	cd->flag = 2;
    }
}

static void fp_mul_si(element_ptr c, element_ptr a, signed long int op)
{
    dataptr ad = a->data;
    dataptr cd = c->data;

    if (!ad->flag || !op) {
	cd->flag = 0;
    } else {
	cd->flag = 2;
	fp_field_data_ptr p = a->field->data;
	size_t t = p->limbs;
	mp_limb_t tmp[t + 1];
	mp_limb_t qp[2];

	tmp[t] = mpn_mul_1(tmp, ad->d, t, labs(op));
	mpn_tdiv_qr(qp, cd->d, 0, tmp, t + 1, p->primelimbs, t);
	if (op < 0) { //TODO: don't need to check c != 0 this time
	    fp_neg(c, c);
	}
    }
}

static void fp_pow_mpz(element_ptr c, element_ptr a, mpz_ptr op)
{
    dataptr ad = a->data;
    dataptr cd = c->data;
    if (!ad->flag) cd->flag = 0;
    else {
	mpz_t z;
	mpz_init(z);
	fp_to_mpz(z, a);
	mpz_powm(z, z, op, a->field->order);
	from_mpz(c, z);
	mpz_clear(z);
	cd->flag = 2;
    }
}

static void fp_invert(element_ptr c, element_ptr a)
{
    //assumes a is invertible
    dataptr cd = c->data;
    mpz_t z;
    mpz_init(z);
    fp_to_mpz(z, a);
    mpz_invert(z, z, a->field->order);
    from_mpz(c, z);
    mpz_clear(z);
    cd->flag = 2;
}

static void fp_random(element_ptr a)
{
    dataptr ad = a->data;
    mpz_t z;
    mpz_init(z);
    pbc_mpz_random(z, a->field->order);
    if (mpz_sgn(z)) {
	from_mpz(a, z);
	ad->flag = 2;
    } else {
	ad->flag = 0;
    }
    mpz_clear(z);
}

static void fp_from_hash(element_ptr a, int len, void *data)
    //TODO: something more sophisticated!
{
    mpz_t z;

    mpz_init(z);
    mpz_import(z, len, -1, 1, -1, 0, data);
    fp_set_mpz(a, z);
    mpz_clear(z);
}

static int fp_cmp(element_ptr a, element_ptr b)
{
    dataptr ad = a->data, bd = b->data;
    if (!ad->flag) {
	return !bd->flag;
    } else {
	fp_field_data_ptr p = a->field->data;
	return mpn_cmp(ad->d, bd->d, p->limbs);
	//return memcmp(ad->d, bd->d, p->limbs);
    }
}

static int fp_is_sqr(element_ptr a)
{
    dataptr ad = a->data;
    int res;
    mpz_t z;
    mpz_init(z);
    //0 is a square
    if (!ad->flag) return 1;
    fp_to_mpz(z, a);
    res = mpz_legendre(z, a->field->order) == 1;
    mpz_clear(z);
    return res;
}

//x is stored as xR, hence apply reduction before output
static int fp_to_bytes(unsigned char *data, element_t a)
{
    dataptr ad = a->data;
    int n = a->field->fixed_length_in_bytes;
    if (!ad->flag) {
	memset(data, 0, n);
    } else {
	fp_field_data_ptr p = a->field->data;
	mp_limb_t tmp[2 * p->limbs];
	mp_limb_t out[p->limbs];

	memcpy(tmp, ad->d, p->limbs * sizeof(mp_limb_t));
	memset(&tmp[p->limbs], 0, p->limbs * sizeof(mp_limb_t));
	mont_reduce(out, tmp, p);
	memcpy(data, out, n);
    }
    return n;
}

static int fp_from_bytes(element_t a, unsigned char *data)
{
    dataptr ad = a->data;
    unsigned char *ptr;
    int i, n;
    mpz_t z, z1;

    mpz_init(z);
    mpz_init(z1);
    mpz_set_ui(z, 0);

    ptr = data;
    n = a->field->fixed_length_in_bytes;
    for (i=0; i<n; i++) {
	if (*ptr) ad->flag = 2;
	mpz_set_ui(z1, *ptr);
	mpz_mul_2exp(z1, z1, i * 8);
	ptr++;
	mpz_add(z, z, z1);
    }
    if (ad->flag) from_mpz(a, z);
    mpz_clear(z);
    mpz_clear(z1);
    return n;
}

static void fp_field_clear(field_t f)
{
    fp_field_data_ptr p = f->data;
    free(p->primelimbs);
    free(p);
}

void field_init_mont_fp(field_ptr f, mpz_t prime)
{
    assert (!mpz_fits_ulong_p(prime));
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
    f->mul_si = fp_mul_si;
    f->doub = fp_double;
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
    {
	mpz_t z;
	mpz_init(z);
	size_t count;

	p->R = malloc(p->bytes);
	mpz_setbit(z, p->bytes * 8);
	mpz_mod(z, z, prime);
	mpz_export(p->R, &count, -1, sizeof(mp_limb_t), 0, 0, z);
	memset((void *) (((unsigned int) p->R) + count * sizeof(mp_limb_t)),
		0, (p->limbs - count) * sizeof(mp_limb_t));
	mpz_set_ui(z, 0);

	//TODO: Algorithm II.5 in Blake, Seroussi and Smart is better
	//but this will do since we're only doing it once
	mpz_setbit(z, p->bytes * 8);
	mpz_invert(z, prime, z);
	p->negpinv = -mpz_get_ui(z);
	mpz_clear(z);
    }
}
