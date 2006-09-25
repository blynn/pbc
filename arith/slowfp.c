//TODO: XXX WARNING XXX this file is very broken
//I can't get mont_reduce faster than mpz_mod() for some reason
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "field.h"
#include "random.h"
//Use Montgomery method for faster multiplication
//note inversion is slower
//an element of F_p is represented by xR (mod p)
//where R is the smallest power of 2 greater than p

//TODO: it's slower than the naive implementation

struct fp_field_data_s {
    int Rwords;
    unsigned long int minuspinv;

    int Rbits;
    mpz_t R;
    mpz_t Rmodp;
    mpz_t R3modp;
};
typedef struct fp_field_data_s fp_field_data_t[1];
typedef struct fp_field_data_s *fp_field_data_ptr;

static void fp_init(element_ptr e)
{
    e->data = malloc(sizeof(mpz_t));
    mpz_init(e->data);
}

static void fp_clear(element_ptr e)
{
    mpz_clear(e->data);
    free(e->data);
}

//Montgomery reduction
//TODO: why is this slower than mpz_mod?
static void mont_reduce(mpz_ptr z, element_ptr y)
{
    fp_field_data_ptr p = y->field->data;
    mpz_ptr q = y->field->order;
    mpz_t z0, y2;
    mpz_init(z0);
    int i;
    //unsigned long int u;

    for (i=0; i<p->Rwords; i++) {
	//u = mpz_get_ui(y2);
	mpz_mul_ui(z0, q, 1); //mpz_getlimbn(y->data, i) * p->minuspinv);
	mpz_mul_2exp(z0, z0, i * sizeof(unsigned long int) * 8);
	mpz_add(y->data, y->data, z0);
    }

    mpz_tdiv_q_2exp(z, y->data, p->Rbits);
    if (mpz_cmp(z, q) > 0) {
	mpz_sub(z, z, q);
    }

    mpz_clear(z0);
    mpz_clear(y2);
}

static void fp_to_mpz(mpz_ptr z, element_ptr y)
{
    mont_reduce(z, y);
}

//replace x with xR (mod P)
static inline void to_mont(element_ptr e)
{
    mpz_ptr z = e->data;
    fp_field_data_ptr p = e->field->data;
    mpz_mul_2exp(z, z, p->Rbits);
    mpz_mod(z, z, e->field->order);
}

static void fp_set_si(element_ptr e, signed long int op)
{
    mpz_set_si(e->data, op);
    to_mont(e);
}

static void fp_set_mpz(element_ptr e, mpz_ptr z)
{
    mpz_set(e->data, z);
    to_mont(e);
}

static void fp_set0(element_ptr e)
{
    mpz_set_si(e->data, 0);
}

static void fp_set1(element_ptr e)
{
    fp_field_data_ptr p = e->field->data;
    mpz_set(e->data, p->Rmodp);
}

static size_t fp_out_str(FILE *stream, int base, element_ptr e)
{
    size_t result;
    mpz_t z;
    mpz_init(z);

    element_to_mpz(z, e);
    result = mpz_out_str(stream, base, z);
    mpz_clear(z);
    return result;
}

static void fp_add(element_ptr n, element_ptr a, element_ptr b)
{
    mpz_add(n->data, a->data, b->data);
    if (mpz_cmp(n->data, n->field->order) >= 0) {
	mpz_sub(n->data, n->data, n->field->order);
    }
}

static void fp_sub(element_ptr n, element_ptr a, element_ptr b)
{
    mpz_sub(n->data, a->data, b->data);
    if (mpz_sgn((mpz_ptr) n->data) < 0) {
	mpz_add(n->data, n->data, n->field->order);
    }
}

static void fp_mul(element_ptr n, element_ptr a, element_ptr b)
{
    mpz_mul(n->data, a->data, b->data);
    mont_reduce(n->data, n);
    mpz_set(n->data, b->data);
}

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
    fp_field_data_ptr p = a->field->data;
    mpz_invert(n->data, a->data, n->field->order);
    mpz_mul(n->data, n->data, p->R3modp);
    mont_reduce(n->data, n);
}

static void fp_random(element_ptr n)
{
    pbc_mpz_random(n->data, n->field->order);
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

static int fp_is1(element_ptr n)
{
    fp_field_data_ptr p = n->field->data;
    return !mpz_cmp((mpz_ptr) n->data, p->Rmodp);
}

static int fp_is0(element_ptr n)
{
    return mpz_is0(n->data);
}

static int fp_cmp(element_ptr a, element_ptr b)
{
    return mpz_cmp((mpz_ptr) a->data, (mpz_ptr) b->data);
}

static int fp_is_sqr(element_ptr a)
{
    mpz_t z;
    int res;

    //0 is a square
    if (mpz_is0(a->data)) return 1;
    mpz_init(z);
    mont_reduce(z, a);
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
    mpz_clear(p->R);
    mpz_clear(p->Rmodp);
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

void field_init_slow_fp(field_ptr f, mpz_t prime)
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
    mpz_init(p->R);
    mpz_init(p->Rmodp);
    mpz_init(p->R3modp);

    int len = sizeof(unsigned long int) * 8;
    mpz_t z;
    mpz_init(z);
    
    p->Rwords = (mpz_sizeinbase(prime, 2) + (len - 1)) / len;
    p->Rbits = p->Rwords * sizeof(unsigned long int) * 8;

    mpz_setbit(z, len);
    mpz_invert(z, prime, z);
    p->minuspinv = -mpz_get_ui(z);

printf("limb: %d\n", sizeof(mp_limb_t));
printf("%d\n",p->Rbits);
mpz_out_str(stdout, 0, prime);
printf("= p\n");
printf("%lu = -p^-1, %lx\n", p->minuspinv, p->minuspinv * mpz_get_ui(prime));
    mpz_mod(p->Rmodp, p->R, prime);
    mpz_mul(p->R3modp, p->Rmodp, p->Rmodp);
    mpz_mod(p->R3modp, p->R3modp, prime);
    mpz_mul(p->R3modp, p->R3modp, p->Rmodp);
    mpz_mod(p->R3modp, p->R3modp, prime);
}
