#include <stdio.h>
#include <stdlib.h>
#include <alloca.h>
#include <string.h>
#include <gmp.h>
#include "field.h"
#include "random.h"
//F_p for small (p at most sizeof(long) bytes long)

static void fp_init(element_ptr e)
{
    e->data = malloc(sizeof(unsigned long));
}

static void fp_clear(element_ptr e)
{
    free(e->data);
}

static void fp_set_mpz(element_ptr e, mpz_ptr z)
{
    mpz_t r;
    mpz_init(r);
    unsigned long *p = e->field->data;
    unsigned long *l = e->data;
    mpz_fdiv_r_ui(r, z, *p);
    *l = mpz_get_ui(r);
    mpz_clear(r);
}

static void fp_set_si(element_ptr e, signed long int op)
{
    unsigned long int *d = e->data;
    unsigned long *p = e->field->data;
    if (op < 0) {
	*d = (-op) % *p;
	*d = *p - *d;
    } else {
	*d = op % *p;
    }
}

static void fp_to_mpz(mpz_ptr z, element_ptr e)
{
    unsigned long int *l = e->data;
    mpz_set_ui(z, *l);
}

static void fp_set0(element_ptr e)
{
    unsigned long int *l = e->data;
    *l = 0;
}

static void fp_set1(element_ptr e)
{
    unsigned long int *l = e->data;
    *l = 1;
}

static int fp_is1(element_ptr e)
{
    unsigned long int *l = e->data;
    return *l == 1;
}

static int fp_is0(element_ptr e)
{
    unsigned long int *l = e->data;
    return *l == 0;
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

static void fp_add(element_ptr c, element_ptr a, element_ptr b)
{
    unsigned long *prime = a->field->data;
    unsigned long *p = a->data;
    unsigned long *q = b->data;
    unsigned long *r = c->data;
    unsigned long l0;
    l0 = *p + *q;
    if (l0 < *p) {
	//overflow
	l0 -= *prime;
    }
    *r = l0 % *prime;
}

static void fp_double(element_ptr c, element_ptr a)
{
    unsigned long *prime = a->field->data;
    unsigned long *p = a->data;
    unsigned long *r = c->data;
    *r = 2 * *p;
    if (*r < *p) {
	//overflow
	*r -= *prime;
    }
    *r = *r % *prime;
}

static void fp_sub(element_ptr c, element_ptr a, element_ptr b)
{
    unsigned long *prime = a->field->data;
    unsigned long *p = a->data;
    unsigned long *q = b->data;
    unsigned long *r = c->data;

    if (*p >= *q) {
	*r = *p - *q;
    } else {
	*r = *prime + *q - *p;
    }
}

static void fp_mul(element_ptr c, element_ptr a, element_ptr b)
{
    unsigned long *prime = a->field->data;
    unsigned long *p = a->data;
    unsigned long *q = b->data;
    unsigned long *r = c->data;

    //TODO: overflow!
    *r = (*p * *q) % *prime;
}

static void fp_square(element_ptr c, element_ptr a)
{
    fp_mul(c, a, a);
}

static void fp_neg(element_ptr c, element_ptr a)
{
    unsigned long *prime = a->field->data;
    unsigned long *r = c->data;
    unsigned long *p = a->data;
    if (*p) {
	*r = *prime - *p;
    } else {
	*r = 0;
    }
}

static void fp_mul_si(element_ptr e, element_ptr a, signed long int op)
{
    //TODO
}

static void fp_pow_mpz(element_ptr n, element_ptr a, mpz_ptr op)
{
    //TODO
    mpz_t z;
    mpz_init(z);
    fp_to_mpz(z, a);
    //mpz_powm(z, z, op, n->field->order);
    //from_mpz(n, z);
    mpz_clear(z);
}

static void fp_set(element_ptr c, element_ptr a)
{
    unsigned long *p = a->data;
    unsigned long *r = c->data;
    *r = *p;
}

static void fp_invert(element_ptr e, element_ptr a)
{
    /* TODO
    mpz_t z;
    mpz_init(z);
    fp_to_mpz(z, a);
    mpz_invert(z, z, e->field->order);
    from_mpz(e, z);
    mpz_clear(z);
    */
}

static void fp_random(element_ptr n)
{
    /* TODO
    mpz_t z;
    mpz_init(z);
    pbc_mpz_random(z, n->field->order);
    from_mpz(n, z);
    mpz_clear(z);
    */
}

static void fp_from_hash(element_ptr n, int len, void *data)
{
    /* TODO
    mpz_t z;

    mpz_init(z);
    mpz_import(z, len, 1, 1, 0, 0, data);
    fp_set_mpz(n, z);
    mpz_clear(z);
    */
}

static int fp_cmp(element_ptr a, element_ptr b)
{
    unsigned long *p = a->data;
    unsigned long *q = b->data;
    return *p != *q;
}

static int fp_is_sqr(element_ptr a)
{
    /* TODO
    int res;
    mpz_t z;
    mpz_init(z);
    //0 is a square
    if (fp_is0(a)) return 1;
    fp_to_mpz(z, a);
    res = mpz_legendre(z, a->field->order) == 1;
    mpz_clear(z);
    return res;
    */
}

static void fp_tonelli(element_ptr x, element_ptr a)
{
    /* TODO
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
    */
    /* TODO: this would be a good place to use element_pow2 ... -hs */
    /*
    element_mul(x, x, e0);
    mpz_clear(t);
    mpz_clear(e);
    mpz_clear(t0);
    element_clear(ginv);
    element_clear(e0);
    */
}

static int fp_to_bytes(unsigned char *data, element_t e)
{
    unsigned char *ptr = data;
    unsigned long *p = e->data;
    unsigned long l = *p;
    int i, n = sizeof(unsigned long);
    for (i = 0; i < n; i++) {
	*ptr = (unsigned char) l;
	ptr++;
	l >>= 8;
    }
    return n;
}

static int fp_from_bytes(element_t e, unsigned char *data)
{
    unsigned char *ptr;
    unsigned long *p = e->data;
    unsigned long l = 0;
    int i, n = sizeof(unsigned long);
    *p = 0;
    for (i=0; i<n; i++) {
	l = *ptr;
	l <<= 8 * i;
	*p += l;
    }
    return n;
}

static void fp_field_clear(field_t f)
{
    free(f->data);
}

void field_init_tiny_fp(field_ptr f, mpz_t prime)
{
    unsigned long *p;

    if (mpz_sizeinbase(prime, 2) > sizeof(long) * 8) {
	printf("field is too big!\n");
	exit(1);
    }
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
    f->square = fp_square;
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

    p = f->data = malloc(sizeof(long));
    p = mpz_get_ui(prime);
    f->fixed_length_in_bytes = sizeof(long);
}
