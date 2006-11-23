#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "pbc_field.h"
#include "pbc_fieldquadratic.h"
#include "pbc_utils.h"
#include "pbc_memory.h"

static inline element_ptr fq_nqr(field_ptr f)
{
    return ((field_ptr) f->data)->nqr;
}

static void fq_init(element_ptr e)
{
    fq_data_ptr p = e->data = pbc_malloc(sizeof(fq_data_t));
    field_ptr f = e->field->data;
    element_init(p->x, f);
    element_init(p->y, f);
}

static void fq_clear(element_ptr e)
{
    fq_data_ptr p = e->data;
    element_clear(p->x);
    element_clear(p->y);
    pbc_free(e->data);
}

static void fq_set_si(element_ptr e, signed long int i)
{
    fq_data_ptr p = e->data;
    element_set_si(p->x, i);
    element_set0(p->y);
}

static void fq_set_mpz(element_ptr e, mpz_t z)
{
    fq_data_ptr p = e->data;
    element_set_mpz(p->x, z);
    element_set0(p->y);
}

//projection: attempts to convert Re(e) to mpz
static void fq_to_mpz(mpz_t z, element_ptr e)
{
    fq_data_ptr p = e->data;
    element_to_mpz(z, p->x);
}

static void fq_set0(element_ptr e)
{
    fq_data_ptr p = e->data;
    element_set0(p->x);
    element_set0(p->y);
}

static void fq_set1(element_ptr e)
{
    fq_data_ptr p = e->data;
    element_set1(p->x);
    element_set0(p->y);
}

static int fq_is0(element_ptr e)
{
    fq_data_ptr p = e->data;
    return element_is0(p->x) && element_is0(p->y);
}

static int fq_is1(element_ptr e)
{
    fq_data_ptr p = e->data;
    return element_is1(p->x) && element_is0(p->y);
}

static size_t fq_out_str(FILE *stream, int base, element_ptr e)
{
    size_t result = 3, status;
    fq_data_ptr p = e->data;
    if (EOF == fputc('[', stream)) return 0;
    result = element_out_str(stream, base, p->x);
    if (!result) return 0;
    if (EOF == fputc(' ', stream)) return 0;
    status = element_out_str(stream, base, p->y);
    if (!status) return 0;
    if (EOF == fputc(']', stream)) return 0;
    return result + status;
}

static int fq_sign(element_ptr n)
{
    int res;
    fq_data_ptr r = n->data;
    res = element_sign(r->x);
    if (!res) return element_sign(r->y);
    return res;
}

static void fq_add(element_ptr n, element_ptr a, element_ptr b)
{
    fq_data_ptr p = a->data;
    fq_data_ptr q = b->data;
    fq_data_ptr r = n->data;
    element_add(r->x, p->x, q->x);
    element_add(r->y, p->y, q->y);
}

static void fq_double(element_ptr n, element_ptr a)
{
    fq_data_ptr p = a->data;
    fq_data_ptr r = n->data;
    element_double(r->x, p->x);
    element_double(r->y, p->y);
}

static void fq_sub(element_ptr n, element_ptr a, element_ptr b)
{
    fq_data_ptr p = a->data;
    fq_data_ptr q = b->data;
    fq_data_ptr r = n->data;
    element_sub(r->x, p->x, q->x);
    element_sub(r->y, p->y, q->y);
}

static void fq_set(element_ptr n, element_ptr a)
{
    fq_data_ptr p = a->data;
    fq_data_ptr r = n->data;
    element_set(r->x, p->x);
    element_set(r->y, p->y);
}

static void fq_mul(element_ptr n, element_ptr a, element_ptr b)
{
    fq_data_ptr p = a->data;
    fq_data_ptr q = b->data;
    fq_data_ptr r = n->data;

    element_ptr nqr = fq_nqr(n->field);
    element_t e0, e1, e2;

    element_init(e0, p->x->field);
    element_init(e1, e0->field);
    element_init(e2, e0->field);
    /* naive:
    element_mul(e0, p->x, q->x);
    element_mul(e1, p->y, q->y);
    element_mul(e1, e1, nqr);
    element_add(e0, e0, e1);
    element_mul(e1, p->x, q->y);
    element_mul(e2, p->y, q->x);
    element_add(e1, e1, e2);
    element_set(r->x, e0);
    element_set(r->y, e1);
    */
    //Karatsuba:
    element_add(e0, p->x, p->y);
    element_add(e1, q->x, q->y);
    element_mul(e2, e0, e1);
    element_mul(e0, p->x, q->x);
    element_mul(e1, p->y, q->y);
    element_mul(r->x, e1, nqr);
    element_add(r->x, r->x, e0);
    element_sub(e2, e2, e0);
    element_sub(r->y, e2, e1);

    element_clear(e0);
    element_clear(e1);
    element_clear(e2);
}

static void fq_mul_mpz(element_ptr n, element_ptr a, mpz_ptr z)
{
    fq_data_ptr p = a->data;
    fq_data_ptr r = n->data;
    element_mul_mpz(r->x, p->x, z);
    element_mul_mpz(r->y, p->y, z);
}

static void fq_mul_si(element_ptr n, element_ptr a, signed long int z)
{
    fq_data_ptr p = a->data;
    fq_data_ptr r = n->data;
    element_mul_si(r->x, p->x, z);
    element_mul_si(r->y, p->y, z);
}

static void fq_square(element_ptr n, element_ptr a)
{
    fq_data_ptr p = a->data;
    fq_data_ptr r = n->data;
    element_ptr nqr = fq_nqr(n->field);
    element_t e0, e1;

    element_init(e0, p->x->field);
    element_init(e1, e0->field);
    element_square(e0, p->x);
    element_square(e1, p->y);
    element_mul(e1, e1, nqr);
    element_add(e0, e0, e1);
    element_mul(e1, p->x, p->y);
    //TODO: which is faster?
    //element_add(e1, e1, e1);
    element_double(e1, e1);
    element_set(r->x, e0);
    element_set(r->y, e1);
    element_clear(e0);
    element_clear(e1);
}

static void fq_neg(element_ptr n, element_ptr a)
{
    fq_data_ptr p = a->data;
    fq_data_ptr r = n->data;
    element_neg(r->x, p->x);
    element_neg(r->y, p->y);
}

static void fq_random(element_ptr e)
{
    fq_data_ptr p = e->data;
    element_random(p->x);
    element_random(p->y);
}

static int fq_cmp(element_ptr a, element_ptr b)
{
    fq_data_ptr p = a->data;
    fq_data_ptr q = b->data;
    return element_cmp(p->x, q->x) || element_cmp(p->y, q->y);
}

static void fq_invert(element_ptr n, element_ptr a)
{
    fq_data_ptr p = a->data;
    fq_data_ptr r = n->data;
    element_ptr nqr = fq_nqr(n->field);
    element_t e0, e1;

    element_init(e0, p->x->field);
    element_init(e1, e0->field);
    element_square(e0, p->x);
    element_square(e1, p->y);
    element_mul(e1, e1, nqr);
    element_sub(e0, e0, e1);
    element_invert(e0, e0);
    element_mul(r->x, p->x, e0);
    element_neg(e0, e0);
    element_mul(r->y, p->y, e0);

    element_clear(e0);
    element_clear(e1);
}

static void fq_from_hash(element_ptr n, void *data, int len)
{
    fq_data_ptr r = n->data;
    int k = len / 2;
    element_from_hash(r->x, data, k);
    element_from_hash(r->y, (char *)data + k, len - k);
}

static int fq_length_in_bytes(element_ptr e)
{
    fq_data_ptr p = e->data;
    return element_length_in_bytes(p->x) + element_length_in_bytes(p->y);
}

static int fq_to_bytes(unsigned char *data, element_t e)
{
    fq_data_ptr p = e->data;
    int len;
    len = element_to_bytes(data, p->x);
    len += element_to_bytes(data + len, p->y);
    return len;
}

static int fq_from_bytes(element_t e, unsigned char *data)
{
    fq_data_ptr p = e->data;
    int len;
    len = element_from_bytes(p->x, data);
    len += element_from_bytes(p->y, data + len);
    return len;
}

static int fq_is_sqr(element_ptr e)
{
    //x + y sqrt(nqr) is a square iff x^2 - nqr y^2 is (in the base field)
    fq_data_ptr p = e->data;
    element_t e0, e1;
    element_ptr nqr = fq_nqr(e->field);
    int result;
    element_init(e0, p->x->field);
    element_init(e1, e0->field);
    element_square(e0, p->x);
    element_square(e1, p->y);
    element_mul(e1, e1, nqr);
    element_sub(e0, e0, e1);
    result = element_is_sqr(e0);
    element_clear(e0);
    element_clear(e1);
    return result;
}

static void fq_sqrt(element_ptr n, element_ptr e)
{
    fq_data_ptr p = e->data;
    fq_data_ptr r = n->data;
    element_ptr nqr = fq_nqr(n->field);
    element_t e0, e1, e2;

    //if (a+b sqrt(nqr))^2 = x+y sqrt(nqr) then
    //2a^2 = x +- sqrt(x^2 - nqr y^2)
    //(take the sign which allows a to exist)
    //and 2ab = y
    element_init(e0, p->x->field);
    element_init(e1, e0->field);
    element_init(e2, e0->field);
    element_square(e0, p->x);
    element_square(e1, p->y);
    element_mul(e1, e1, nqr);
    element_sub(e0, e0, e1);
    element_sqrt(e0, e0);
    //e0 = sqrt(x^2 - nqr y^2)
    element_add(e1, p->x, e0);
    element_set_si(e2, 2);
    element_invert(e2, e2);
    element_mul(e1, e1, e2);
    //e1 = (x + sqrt(x^2 - nqr y^2))/2
    if (!element_is_sqr(e1)) {
	element_sub(e1, e1, e0);
	//e1 should be a square
    }
    element_sqrt(e0, e1);
    element_add(e1, e0, e0);
    element_invert(e1, e1);
    element_mul(r->y, p->y, e1);
    element_set(r->x, e0);
    element_clear(e0);
    element_clear(e1);
    element_clear(e2);
}

static void field_clear_fq(field_ptr f)
{
    UNUSED_VAR(f);
    //f->order gets cleared automatically
}

static void fq_out_info(FILE *out, field_ptr f)
{
    field_ptr fbase = f->data;
    element_fprintf(out, "x^2 + %B quadratic extension, base field:\n", fq_nqr(f));
    field_out_info(out, fbase);
}

void field_init_quadratic(field_ptr f, field_ptr fbase)
{
    field_init(f);

    f->field_clear = field_clear_fq;
    f->data = fbase;

    f->init = fq_init;
    f->clear = fq_clear;
    f->set_si = fq_set_si;
    f->set_mpz = fq_set_mpz;
    f->to_mpz = fq_to_mpz;
    f->out_str = fq_out_str;
    f->sign = fq_sign;
    f->add = fq_add;
    f->sub = fq_sub;
    f->set = fq_set;
    f->mul = fq_mul;
    f->mul_mpz = fq_mul_mpz;
    f->mul_si = fq_mul_si;
    f->square = fq_square;
    f->doub = fq_double;
    f->neg = fq_neg;
    f->cmp = fq_cmp;
    f->invert = fq_invert;
    f->random = fq_random;
    f->from_hash = fq_from_hash;
    f->is1 = fq_is1;
    f->is0 = fq_is0;
    f->set0 = fq_set0;
    f->set1 = fq_set1;
    f->is_sqr = fq_is_sqr;
    f->sqrt = fq_sqrt;
    f->to_bytes = fq_to_bytes;
    f->from_bytes = fq_from_bytes;
    f->out_info = fq_out_info;

    mpz_mul(f->order, fbase->order, fbase->order);
    if (fbase->fixed_length_in_bytes < 0) {
	f->length_in_bytes = fq_length_in_bytes;
	f->fixed_length_in_bytes = -1;
    } else {
	f->fixed_length_in_bytes = 2 * fbase->fixed_length_in_bytes;
    }
}

void element_field_to_quadratic(element_ptr r, element_ptr a)
{
    fq_data_ptr p = r->data;
    element_set(p->x, a);
    element_set0(p->y);
}

static void fi_mul(element_ptr n, element_ptr a, element_ptr b)
{
    fq_data_ptr p = a->data;
    fq_data_ptr q = b->data;
    fq_data_ptr r = n->data;
    element_t e0, e1, e2;

    element_init(e0, p->x->field);
    element_init(e1, e0->field);
    element_init(e2, e0->field);
    /* Naive way
    element_mul(e0, p->x, q->x);
    element_mul(e1, p->y, q->y);
    element_sub(e0, e0, e1);
    element_mul(e1, p->x, q->y);
    element_mul(e2, p->y, q->x);
    element_add(e1, e1, e2);
    element_set(r->x, e0);
    element_set(r->y, e1);
    */
    //Karatsuba:
    element_add(e0, p->x, p->y);
    element_add(e1, q->x, q->y);
    element_mul(e2, e0, e1);
    element_mul(e0, p->x, q->x);
    element_sub(e2, e2, e0);
    element_mul(e1, p->y, q->y);
    element_sub(r->x, e0, e1);
    element_sub(r->y, e2, e1);

    element_clear(e0);
    element_clear(e1);
    element_clear(e2);
}

static void fi_square(element_ptr n, element_ptr a)
{
    fq_data_ptr p = a->data;
    fq_data_ptr r = n->data;
    element_t e0, e1;

    element_init(e0, p->x->field);
    element_init(e1, e0->field);
    //Re(n) = x^2 - y^2 = (x+y)(x-y)
    element_add(e0, p->x, p->y);
    element_sub(e1, p->x, p->y);
    element_mul(e0, e0, e1);
    //Im(n) = 2xy
    element_mul(e1, p->x, p->y);
    element_add(e1, e1, e1);
    element_set(r->x, e0);
    element_set(r->y, e1);
    element_clear(e0);
    element_clear(e1);
}

static void fi_invert(element_ptr n, element_ptr a)
{
    fq_data_ptr p = a->data;
    fq_data_ptr r = n->data;
    element_t e0, e1;

    element_init(e0, p->x->field);
    element_init(e1, e0->field);
    element_square(e0, p->x);
    element_square(e1, p->y);
    element_add(e0, e0, e1);
    element_invert(e0, e0);
    element_mul(r->x, p->x, e0);
    element_neg(e0, e0);
    element_mul(r->y, p->y, e0);

    element_clear(e0);
    element_clear(e1);
}

static int fi_is_sqr(element_ptr e)
{
    //x + yi is a square <=> x^2 + y^2 is (in the base field)

    // Proof: (=>) if x+yi = (a+bi)^2,
    // then a^2 - b^2 = x, 2ab = y,
    // thus (a^2 + b^2)^2 = (a^2 - b^2)^2 + (2ab)^2 =  x^2 + y^2
    // (<=) Suppose A^2 = x^2 + y^2
    // then if there exist a, b satisfying:
    //   a^2 = (+-A + x)/2, b^2 = (+-A - x)/2
    // then (a + bi)^2 = x + yi.
    // We show that exactly one of (A + x)/2, (-A + x)/2
    // is a quadratic residue (thus a, b do exist).
    // Suppose not. Then the product
    // (x^2 - A^2) / 4 is some quadratic residue, a contradiction
    // since this would imply x^2 - A^2 = -y^2 is also a quadratic residue,
    // but we know -1 is not a quadratic residue.
    fq_data_ptr p = e->data;
    element_t e0, e1;
    int result;
    element_init(e0, p->x->field);
    element_init(e1, e0->field);
    element_square(e0, p->x);
    element_square(e1, p->y);
    element_add(e0, e0, e1);
    result = element_is_sqr(e0);
    element_clear(e0);
    element_clear(e1);
    return result;
}

static void fi_sqrt(element_ptr n, element_ptr e)
{
    fq_data_ptr p = e->data;
    fq_data_ptr r = n->data;
    element_t e0, e1, e2;

    //if (a+bi)^2 = x+yi then
    //2a^2 = x +- sqrt(x^2 + y^2)
    //(take the sign such that a exists) and 2ab = y
    //[thus 2b^2 = - (x -+ sqrt(x^2 + y^2))]
    element_init(e0, p->x->field);
    element_init(e1, e0->field);
    element_init(e2, e0->field);
    element_square(e0, p->x);
    element_square(e1, p->y);
    element_add(e0, e0, e1);
    element_sqrt(e0, e0);
    //e0 = sqrt(x^2 + y^2)
    element_add(e1, p->x, e0);
    element_set_si(e2, 2);
    element_invert(e2, e2);
    element_mul(e1, e1, e2);
    //e1 = (x + sqrt(x^2 + y^2))/2
    if (!element_is_sqr(e1)) {
	element_sub(e1, e1, e0);
	//e1 should be a square
    }
    element_sqrt(e0, e1);
    element_add(e1, e0, e0);
    element_invert(e1, e1);
    element_mul(r->y, p->y, e1);
    element_set(r->x, e0);
    element_clear(e0);
    element_clear(e1);
    element_clear(e2);
}

void element_field_to_fi(element_ptr a, element_ptr b)
{
    fq_data_ptr p = a->data;

    element_set(p->x, b);
    element_set0(p->y);
}

static void fi_out_info(FILE *out, field_ptr f)
{
    field_ptr fbase = f->data;
    fprintf(out, "x^2 + 1 quadratic extension, base field:\n");
    field_out_info(out, fbase);
}

static void field_clear_fi(field_ptr f)
{
    UNUSED_VAR(f);
    //f->order gets cleared automatically
}

void field_init_fi(field_ptr f, field_ptr fbase)
{
    field_init(f);
    f->field_clear = field_clear_fi;
    f->data = fbase;
    f->init = fq_init;
    f->clear = fq_clear;
    f->set_si = fq_set_si;
    f->set_mpz = fq_set_mpz;
    f->to_mpz = fq_to_mpz;
    f->out_str = fq_out_str;
    f->sign = fq_sign;
    f->add = fq_add;
    f->sub = fq_sub;
    f->set = fq_set;
    f->mul = fi_mul;
    f->mul_mpz = fq_mul_mpz;
    f->mul_si = fq_mul_si;
    f->square = fi_square;
    f->doub = fq_double;
    f->neg = fq_neg;
    f->cmp = fq_cmp;
    f->invert = fi_invert;
    f->random = fq_random;
    f->from_hash = fq_from_hash;
    f->is1 = fq_is1;
    f->is0 = fq_is0;
    f->set0 = fq_set0;
    f->set1 = fq_set1;
    f->is_sqr = fi_is_sqr;
    f->sqrt = fi_sqrt;
    f->to_bytes = fq_to_bytes;
    f->from_bytes = fq_from_bytes;
    f->out_info = fi_out_info;

    mpz_mul(f->order, fbase->order, fbase->order);
    if (fbase->fixed_length_in_bytes < 0) {
	f->length_in_bytes = fq_length_in_bytes;
	f->fixed_length_in_bytes = -1;
    } else {
	f->fixed_length_in_bytes = 2 * fbase->fixed_length_in_bytes;
    }
}
