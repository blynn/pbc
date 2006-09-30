#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "field.h"

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

static void generic_print_info(FILE *out, field_ptr f)
{
    element_fprintf(out, "field %X unknown\n", f);
}

void field_print_info(FILE *out, field_ptr f)
{
    f->print_info(out, f);
}

void field_init(field_ptr f)
{
    f->nqr = NULL;
    mpz_init(f->order);
    f->square = generic_square;
    f->mul_mpz = generic_mul_mpz;
    f->pow_mpz = generic_pow_mpz;
    f->mul_si = generic_mul_si;
    f->print_info = generic_print_info;
    f->field_clear = NULL;
}

void field_clear(field_ptr f)
{
    if (f->nqr) {
	element_clear(f->nqr);
	free(f->nqr);
    }
    mpz_clear(f->order);
    if (f->field_clear) f->field_clear(f);
}
