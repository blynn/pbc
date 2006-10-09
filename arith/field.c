#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> //for memcmp()
#include <gmp.h>
#include "field.h"
#include "utils.h"

/* returns recommended window size.  n is exponent. */
static int optimal_pow_window_size(mpz_ptr n)
{
  int exp_bits;

  exp_bits = mpz_sizeinbase(n, 2);

  /* try to minimize 2^k + n/(k+1). */
  if (exp_bits > 9065)
      return 8;
  if (exp_bits > 3529)
      return 7;
  if (exp_bits > 1324)
      return 6;
  if (exp_bits > 474)
      return 5;
  if (exp_bits > 157)
      return 4;
  if (exp_bits > 47)
      return 3;
  return 2;
}  

/* builds k-bit lookup window for base a */
static element_t *build_pow_window(element_ptr a, int k)
{
    int s;
    int lookup_size;
    element_t *lookup;

    if (k < 1) {                /* no window */
        return NULL;
    }

    /* build 2^k lookup table.  lookup[i] = x^i. */
    /* TODO: a more careful word-finding algorithm would allow
     *       us to avoid calculating even lookup entries > 2
     */
    lookup_size = 1 << k;
    lookup = malloc(lookup_size * sizeof(element_t));
    element_init(lookup[0], a->field);
    element_set1(lookup[0]);
    for (s = 1; s < lookup_size; s++) {
        element_init(lookup[s], a->field);
        element_mul(lookup[s], lookup[s-1], a);
    }

    return lookup;
}

static void clear_pow_window(int k, element_t *lookup)
{
    int s;
    int lookup_size = 1 << k;

    for (s = 0; s < lookup_size; s++) {
        element_clear(lookup[s]);
    }
    free(lookup);
}

/*
 * left-to-right exponentiation with k-bit window.
 * NB. must have k >= 1.
 */
static void element_pow_wind(element_ptr x, mpz_ptr n,
                             int k, element_t *a_lookup)
{
    int s;
    int bit;

    int inword;                 /* boolean: currently reading word? */
    int word = 0;               /* the word to look up. 0<word<base */
    int wbits = 0;           /* # of bits so far in word. wbits<=k. */

    element_t result;


    /* early abort if raising to power 0 */
    if (!mpz_sgn(n)) {
        element_set1(x);
        return;
    }

    element_init(result, x->field);
    element_set1(result);

    for (inword = 0, s = mpz_sizeinbase(n, 2) - 1; s >=0; s--) {
        element_square(result, result);
        bit = mpz_tstbit(n, s);

        if (!inword && !bit)
            continue;           /* keep scanning.  note continue. */

        if (!inword) {          /* was scanning, just found word */
            inword = 1;         /* so, start new word */
            word = 1; wbits = 1;
        } else {
            word = (word << 1) + bit; wbits++; /* continue word */
        }

        if (wbits == k || s == 0) {
            element_mul(result, result, a_lookup[word]);
            inword = 0;
        }
    }

    element_set(x, result);
    element_clear(result);
}

static void generic_pow_mpz(element_ptr x, element_ptr a, mpz_ptr n)
{
    int k;
    element_t *a_lookup;

    if (mpz_is0(n)) {
        element_set1(x);
        return;
    }

    k = optimal_pow_window_size(n);
    a_lookup = build_pow_window(a, k);
    element_pow_wind(x, n, k, a_lookup);
    clear_pow_window(k, a_lookup);
}

void naive_generic_pow_mpz(element_ptr x, element_ptr a, mpz_ptr n)
{
    int s;

    element_t result;

    if (mpz_is0(n)) {
        element_set1(x);
        return;
    }

    element_init(result, x->field);
    element_set1(result);

    for (s = mpz_sizeinbase(n, 2) - 1; s >=0; s--) {
	element_square(result, result);
	if (mpz_tstbit(n, s)) {
	    element_mul(result, result, a);
	}
    }
    element_set(x, result);
    element_clear(result);
}

void element_pow2_mpz(element_ptr x, element_ptr a1, mpz_ptr n1,
                                 element_ptr a2, mpz_ptr n2)
{
    int s, s1, s2;
    int b1, b2;

    element_t result, a1a2;

    if (mpz_is0(n1) && mpz_is0(n2)) {
        element_set1(x);
        return;
    }

    element_init(result, x->field);
    element_set1(result);

    element_init(a1a2, x->field);
    element_mul(a1a2, a1, a2);

    s1 = mpz_sizeinbase(n1, 2) - 1;
    s2 = mpz_sizeinbase(n2, 2) - 1;
    for (s = (s1 > s2) ? s1 : s2; s >=0; s--) {
        element_mul(result, result, result);
        b1 = mpz_tstbit(n1, s); b2 = mpz_tstbit(n2, s);
        if (b1 && b2) {
            element_mul(result, result, a1a2);
        } else if (b1) {
            element_mul(result, result, a1);
        } else if (b2) {
            element_mul(result, result, a2);
        }
    }

    element_set(x, result);
    element_clear(result);
    element_clear(a1a2);
}

void element_pow3_mpz(element_ptr x, element_ptr a1, mpz_ptr n1,
                                 element_ptr a2, mpz_ptr n2,
                                 element_ptr a3, mpz_ptr n3)
{
    int s, s1, s2, s3;
    int b;
    int i;

    element_t result;
    element_t lookup[8];

    if (mpz_is0(n1) && mpz_is0(n2) && mpz_is0(n3)) {
        element_set1(x);
        return;
    }

    element_init(result, x->field);
    element_set1(result);

    for (i=0; i<8; i++)
        element_init(lookup[i], x->field);

    /* build lookup table. */
    element_set1(lookup[0]);
    element_set(lookup[1], a1);
    element_set(lookup[2], a2);
    element_set(lookup[4], a3);
    element_mul(lookup[3], a1, a2);
    element_mul(lookup[5], a1, a3);
    element_mul(lookup[6], a2, a3);
    element_mul(lookup[7], lookup[6], a1);

    /* calculate largest exponent bitsize */
    s1 = mpz_sizeinbase(n1, 2) - 1;
    s2 = mpz_sizeinbase(n2, 2) - 1;
    s3 = mpz_sizeinbase(n3, 2) - 1;
    s = (s1 > s2) ? ((s1 > s3) ? s1 : s3)
                    : ((s2 > s3) ? s2 : s3);

    for (; s >=0; s--) {
        element_mul(result, result, result);
        b = (mpz_tstbit(n1, s))
          + (mpz_tstbit(n2, s) << 1)
          + (mpz_tstbit(n3, s) << 2);
        element_mul(result, result, lookup[b]);
    }
    
    element_set(x, result);
    element_clear(result);
    for (i=0; i<8; i++)
        element_clear(lookup[i]);
}

struct element_base_table {
    int k;
    int bits;
    element_t **table;
};

/* build k-bit base table for n-bit exponentiation w/ base a */
static void *element_build_base_table(element_ptr a, int bits, int k)
{
    struct element_base_table *base_table;
    element_t multiplier;
    int i, j;
    int lookup_size;
    int num_lookups;

    element_t *lookup;

    fprintf(stderr, "building %d bits %d k\n", bits, k);

    lookup_size = 1 << k;
    num_lookups = bits/k + 1;

    base_table = malloc(sizeof(struct element_base_table));
    base_table->k = k;
    base_table->bits = bits;
    base_table->table = malloc(num_lookups * sizeof(element_t *));

    element_init(multiplier, a->field);
    element_set(multiplier, a);

    for (i = 0; i < num_lookups; i++) {
        lookup = malloc(lookup_size * sizeof(element_t));
        element_init(lookup[0], a->field);
        element_set1(lookup[0]);
        for (j = 1; j < lookup_size; j++) {
            element_init(lookup[j], a->field);
            element_mul(lookup[j], multiplier, lookup[j-1]);
        }
        element_mul(multiplier, multiplier, lookup[lookup_size-1]);
        base_table->table[i] = lookup;
    }

    element_clear(multiplier);
    return base_table;
}

/*
 * exponentiation using aggressive base lookup table
 * must have k >= 1.
 */
static void element_pow_base_table(element_ptr x, mpz_ptr n,
                       struct element_base_table *base_table)
{
    int word;                   /* the word to look up. 0<word<base */
    int row, s;                 /* row and col in base table */
    int num_lookups;

    element_t result;

    /* early abort if raising to power 0 */
    if (!mpz_sgn(n)) {
        element_set1(x);
        return;
    }

    element_init(result, x->field);
    element_set1(result);

    num_lookups = mpz_sizeinbase(n, 2)/base_table->k + 1;

    for (row = 0; row < num_lookups; row++) {
        word = 0;
        for (s = 0; s < base_table->k; s++) {
          word |= mpz_tstbit(n, base_table->k * row + s) << s;
        }
        if (word > 0) {
            element_mul(result, result, base_table->table[row][word]);
        }
    }

    element_set(x, result);
    element_clear(result);
}

void default_element_pp_init(element_pp_t p, element_t in) {
    p->data =
       element_build_base_table(in, mpz_sizeinbase(in->field->order, 2), 5);
}

void default_element_pp_pow(element_t out, mpz_ptr power, element_pp_t p)
{
    element_pow_base_table(out, power, p->data);
}

void default_element_pp_clear(element_pp_t p)
{
    //TODO: free space taken by p->data
    UNUSED_VAR(p);
}

void field_set_nqr(field_ptr f, element_t nqr)
{
    if (!f->nqr) {
	f->nqr = malloc(sizeof(element_t));
	element_init(f->nqr, f);
    }
    element_set(f->nqr, nqr);
}

void field_gen_nqr(field_ptr f)
{
    f->nqr = malloc(sizeof(element_t));
    element_init(f->nqr, f);
    do {
	element_random(f->nqr);
    } while (element_is_sqr(f->nqr));
}

element_ptr field_get_nqr(field_ptr f)
{
    if (!f->nqr) field_gen_nqr(f);
    return f->nqr;
}

static void generic_square(element_ptr r, element_ptr a)
{
    element_mul(r, a, a);
}
static void generic_mul_mpz(element_ptr r, element_ptr a, mpz_ptr z)
{
    element_t e0;
    element_init(e0, r->field);
    element_set_mpz(e0, z);
    element_mul(r, a, e0);
    element_clear(e0);
}

static void generic_mul_si(element_ptr r, element_ptr a, signed long int n)
{
    element_t e0;
    element_init(e0, r->field);
    element_set_si(e0, n);
    element_mul(r, a, e0);
    element_clear(e0);
}

static void generic_double(element_ptr r, element_ptr a)
{
    element_add(r, a, a);
}

static void generic_halve(element_ptr r, element_ptr a)
{
    element_t e0;
    element_init(e0, r->field);
    element_set_si(e0, 2);
    element_invert(e0, e0);
    element_mul(r, a, e0);
    element_clear(e0);
}

static void zero_to_mpz(mpz_t z, element_ptr a)
{
    UNUSED_VAR(a);
    mpz_set_ui(z, 0);
}

static void zero_set_mpz(element_ptr a, mpz_t z)
{
    UNUSED_VAR(z);
    element_set0(a);
}

static void zero_random(element_ptr a)
{
    element_set0(a);
}

static void generic_set_si(element_ptr a, long int si)
{
    mpz_t z;
    mpz_init(z);
    mpz_set_si(z, si);
    element_set_mpz(a, z);
    mpz_clear(z);
}

static void generic_sub(element_ptr c, element_ptr a, element_ptr b)
{
    if (c != a) {
	element_neg(c, b);
	element_add(c, c, a);
    } else {
	element_t tmp;
	element_init(tmp, a->field);
	element_neg(tmp, b);
	element_add(c, tmp, a);
	element_clear(tmp);
    }
}

static void generic_add_ui(element_ptr c, element_ptr a, unsigned long int b)
{
    element_t e;
    mpz_t z;
    element_init(e, c->field);
    mpz_init(z);
    mpz_set_ui(z, b);
    element_set_mpz(e, z);
    element_add(c, a, e);
    mpz_clear(z);
    element_clear(e);
}

static int generic_cmp(element_ptr a, element_ptr b)
{
    int result;
    unsigned char *buf1, *buf2;
    int len;
    if (a == b) return 0;
    len = element_length_in_bytes(a);
    if (len != element_length_in_bytes(b)) return 1;
    buf1 = malloc(len);
    buf2 = malloc(len);
    element_to_bytes(buf1, a);
    element_to_bytes(buf2, b);
    result = memcmp(buf1, buf2, len);
    free(buf1);
    free(buf2);
    return result;
}

static int generic_is0(element_ptr a)
{
    int result;
    element_t b;
    element_init(b, a->field);
    result = element_cmp(a, b);
    element_clear(b);
    return result;
}

static int generic_is1(element_ptr a)
{
    int result;
    element_t b;
    element_init(b, a->field);
    element_set1(b);
    result = element_cmp(a, b);
    element_clear(b);
    return result;
}

static void generic_out_info(FILE *out, field_ptr f)
{
    element_fprintf(out, "field %X unknown\n", (unsigned int) f);
    element_fprintf(out, "order = %Zd\n", f->order);
}

static void warn_field_clear(field_ptr f)
{
    fprintf(stderr, "field %X has no clear function\n", (unsigned int) f);
}

void field_out_info(FILE *out, field_ptr f)
{
    f->out_info(out, f);
}

void field_init(field_ptr f)
{
    //should be called by each field_init_*
    f->nqr = NULL;
    mpz_init(f->order);

    //this should later be set
    f->field_clear = warn_field_clear;

    //and this to something more helpful
    f->out_info = generic_out_info;

    //many of these can usually be optimized for particular fields
    //provided for developer's convenience
    f->halve = generic_halve;
    f->doub = generic_double;
    f->square = generic_square;
    f->mul_mpz = generic_mul_mpz;
    f->mul_si = generic_mul_si;
    f->cmp = generic_cmp;
    f->sub = generic_sub;
    f->add_ui = generic_add_ui;

    //default: converts all elements to integer 0
    //reads all integers as 0
    //random always outputs 0
    f->to_mpz = zero_to_mpz;
    f->set_mpz = zero_set_mpz;
    f->random = zero_random;
    f->set_si = generic_set_si;
    f->is1 = generic_is1;
    f->is0 = generic_is0;

    //these are fast, thanks to Hovav
    f->pow_mpz = generic_pow_mpz;

    f->pp_init = default_element_pp_init;
    f->pp_clear = default_element_pp_clear;
    f->pp_pow = default_element_pp_pow;
}

void field_clear(field_ptr f)
{
    if (f->nqr) {
	element_clear(f->nqr);
	free(f->nqr);
    }
    mpz_clear(f->order);
    f->field_clear(f);
}

void pbc_mpz_out_raw_n(unsigned char *data, int n, mpz_t z)
{
    size_t count;
    if (mpz_sgn(z)) {
	count = (mpz_sizeinbase(z, 2) + 7) / 8;
	mpz_export(&data[n - count], NULL, 1, 1, 1, 0, z);
	memset(data, 0, n - count);
    } else {
	memset(data, 0, n);
    }
}

void pbc_mpz_from_hash(mpz_t z, mpz_t limit,
	unsigned char *data, unsigned int len)
{
    size_t i = 0, n, count = (mpz_sizeinbase(limit, 2) + 7) / 8;
    unsigned char buf[count];
    unsigned char counter = 0;
    int done = 0;
    for(;;) {
	if (len >= count - i) {
	    n = count - i;
	    done = 1;
	} else n = len;
	memcpy(buf + i, data, n);
	i += n;
	if (done) break;
	buf[i] = counter;
	counter++;
	i++;
	if (i == count) break;
    }
    assert(i == count);
    mpz_import(z, count, 1, 1, 1, 0, buf);
    while (mpz_cmp(z, limit) > 0) {
	mpz_tdiv_q_2exp(z, z, 1);
    }
}
