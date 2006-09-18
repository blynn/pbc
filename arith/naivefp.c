#include "field.h"
#include "utils.h"
//Naive implementation of F_p
//may be preferable in some situations

static void zp_init(element_ptr e)
{
    e->data = malloc(sizeof(mpz_t));
    mpz_init(e->data);
}

static void zp_clear(element_ptr e)
{
    mpz_clear(e->data);
    free(e->data);
}

static void zp_set_si(element_ptr e, signed long int op)
{
    mpz_set_si(e->data, op);
    mpz_mod(e->data, e->data, e->field->order);
}

static void zp_set_mpz(element_ptr e, mpz_ptr z)
{
    mpz_set(e->data, z);
    mpz_mod(e->data, e->data, e->field->order);
}

static void zp_set0(element_ptr e)
{
    mpz_set_si(e->data, 0);
}

static void zp_set1(element_ptr e)
{
    mpz_set_si(e->data, 1);
}

static size_t zp_out_str(FILE *stream, int base, element_ptr e)
{
    return mpz_out_str(stream, base, e->data);
}

static int zp_sign(element_ptr a)
{
    mpz_t z;
    mpz_init(z);
    int res;

    mpz_add(z, a->data, a->data);
    if (mpz_is0(z)) {
	res = 0;
    } else {
	res = mpz_cmp(z, a->field->order);
    }
    mpz_clear(z);
    return res;
}

static void zp_add(element_ptr n, element_ptr a, element_ptr b)
{
    /*
    mpz_add(n->data, a->data, b->data);
    mpz_mod(n->data, n->data, n->field->order);
    */
    //This seems faster:
    mpz_add(n->data, a->data, b->data);
    if (mpz_cmp(n->data, n->field->order) >= 0) {
	mpz_sub(n->data, n->data, n->field->order);
    }
}

static void zp_sub(element_ptr n, element_ptr a, element_ptr b)
{
    //mpz_sub(n->data, a->data, b->data);
    //mpz_mod(n->data, n->data, n->field->order);
    mpz_sub(n->data, a->data, b->data);
    if (mpz_sgn((mpz_ptr) n->data) < 0) {
	mpz_add(n->data, n->data, n->field->order);
    }
}

static void zp_square(element_ptr n, element_ptr a)
{
    mpz_mul(n->data, a->data, a->data);
    mpz_mod(n->data, n->data, n->field->order);
}

static void zp_mul(element_ptr n, element_ptr a, element_ptr b)
{
    mpz_mul(n->data, a->data, b->data);
    mpz_mod(n->data, n->data, n->field->order);
}

static void zp_mul_mpz(element_ptr n, element_ptr a, mpz_ptr z)
{
    mpz_mul(n->data, a->data, z);
    mpz_mod(n->data, n->data, n->field->order);
}

static void zp_mul_si(element_ptr n, element_ptr a, signed long int z)
{
    mpz_mul_si(n->data, a->data, z);
    mpz_mod(n->data, n->data, n->field->order);
}

static void zp_pow(element_ptr n, element_ptr a, mpz_ptr z)
{
    mpz_powm(n->data, a->data, z, n->field->order);
}

static void zp_set(element_ptr n, element_ptr a)
{
    mpz_set(n->data, a->data);
}

static void zp_neg(element_ptr n, element_ptr a)
{
    if (mpz_is0(a->data)) {
	mpz_set_ui(n->data, 0);
    } else {
	mpz_sub(n->data, n->field->order, a->data);
    }
}

static void zp_invert(element_ptr n, element_ptr a)
{
    mpz_invert(n->data, a->data, n->field->order);
}

static void zp_random(element_ptr n)
{
    pbc_mpz_random(n->data, n->field->order);
}

static void zp_from_hash(element_ptr n, int len, void *data)
    //TODO: something more sophisticated!
{
    mpz_t z;

    mpz_init(z);
    mpz_import(z, len, 1, 1, 0, 0, data);
    zp_set_mpz(n, z);
    mpz_clear(z);
}

static int zp_is1(element_ptr n)
{
    return !mpz_cmp_ui((mpz_ptr) n->data, 1);
}

static int zp_is0(element_ptr n)
{
    return mpz_is0(n->data);
}

static int zp_cmp(element_ptr a, element_ptr b)
{
    return mpz_cmp((mpz_ptr) a->data, (mpz_ptr) b->data);
}

static int zp_is_sqr(element_ptr a)
{
    //0 is a square
    if (mpz_is0(a->data)) return 1;
    return mpz_legendre(a->data, a->field->order) == 1;
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
	element_pow(e0, ginv, e);
	element_mul(e0, e0, a);
	element_pow(e0, e0, t0);
	if (!element_is1(e0)) mpz_setbit(e, i-1);
    }
    element_pow(e0, ginv, e);
    element_mul(e0, e0, a);
    mpz_add_ui(t, t, 1);
    mpz_tdiv_q_2exp(t, t, 1);
    element_pow(e0, e0, t);
    mpz_tdiv_q_2exp(e, e, 1);
    element_pow(x, nqr, e);
    /* TODO: this would be a good place to use element_pow2 ... -hs */
    element_mul(x, x, e0);
    mpz_clear(t);
    mpz_clear(e);
    mpz_clear(t0);
    element_clear(ginv);
    element_clear(e0);
}

static void zp_field_clear(field_t f)
{
    UNUSED_VAR (f);
}

static int zp_to_bytes(unsigned char *data, element_t e)
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

static int zp_from_bytes(element_t e, unsigned char *data)
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

static void zp_to_mpz(mpz_ptr z, element_ptr a)
{
    mpz_set(z, a->data);
}

void field_init_naive_fp(field_ptr f, mpz_t prime)
{
    field_init(f);
    f->init = zp_init;
    f->clear = zp_clear;
    f->set_si = zp_set_si;
    f->set_mpz = zp_set_mpz;
    f->out_str = zp_out_str;
    f->sign = zp_sign;
    f->add = zp_add;
    f->sub = zp_sub;
    f->set = zp_set;
    f->square = zp_square;
    f->mul = zp_mul;
    f->mul_mpz = zp_mul_mpz;
    f->mul_si = zp_mul_si;
    f->pow = zp_pow;
    f->neg = zp_neg;
    f->cmp = zp_cmp;
    f->invert = zp_invert;
    f->random = zp_random;
    f->from_hash = zp_from_hash;
    f->is1 = zp_is1;
    f->is0 = zp_is0;
    f->set0 = zp_set0;
    f->set1 = zp_set1;
    f->is_sqr = zp_is_sqr;
    f->sqrt = fp_tonelli;
    f->field_clear = zp_field_clear;
    f->to_bytes = zp_to_bytes;
    f->from_bytes = zp_from_bytes;
    f->to_mpz = zp_to_mpz;

    mpz_set(f->order, prime);
    f->data = NULL;
    f->fixed_length_in_bytes = (mpz_sizeinbase(prime, 2) + 7) / 8;
}
