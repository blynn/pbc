// Brute force and Pollard rho discrete log algorithms.

#include <stdarg.h>
#include <stdint.h> // for intptr_t
#include <stdio.h>
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_field.h"
#include "pbc_memory.h"
#include "misc/darray.h"

struct snapshot_s {
  element_t a;
  element_t b;
  element_t snark;
};
typedef struct snapshot_s *snapshot_ptr;

static void record(element_t asum, element_t bsum, element_t snark,
                   darray_t hole, mpz_t counter) {
  snapshot_ptr ss = pbc_malloc(sizeof(struct snapshot_s));
  element_init_same_as(ss->a, asum);
  element_init_same_as(ss->b, bsum);
  element_init_same_as(ss->snark, snark);
  element_set(ss->a, asum);
  element_set(ss->b, bsum);
  element_set(ss->snark, snark);
  darray_append(hole, ss);
  element_printf("snark %Zd: %B\n", counter, snark);
}

// g, h in some group of order r
// finds x such that g^x = h
// will hang if no such x exists
// x in some field_t that set_mpz makes sense for
void element_dlog_brute_force(element_t x, element_t g, element_t h) {
  element_t g0;
  mpz_t count;

  mpz_init(count);
  element_init_same_as(g0, g);

  element_set(g0, g);
  mpz_set_ui(count, 1);
  while (element_cmp(g0, h)) {
    element_mul(g0, g0, g);
//element_printf("g0^%Zd = %B\n", count, g0);
    mpz_add_ui(count, count, 1);
  }
  element_set_mpz(x, count);
  mpz_clear(count);
  element_clear(g0);
}

// x in Z_r, g, h in some group of order r
// finds x such that g^x = h
void element_dlog_pollard_rho(element_t x, element_t g, element_t h) {
// see Blake, Seroussi and Smart
// only one snark for this implementation
  int i, s = 20;
  field_ptr Zr = x->field, G = g->field;
  element_t asum;
  element_t bsum;
  element_t a[s];
  element_t b[s];
  element_t m[s];
  element_t g0, snark;
  darray_t hole;
  int interval = 5;
  mpz_t counter;
  int found = 0;

  mpz_init(counter);
  element_init(g0, G);
  element_init(snark, G);
  element_init(asum, Zr);
  element_init(bsum, Zr);
  darray_init(hole);
  //set up multipliers
  for (i = 0; i < s; i++) {
    element_init(a[i], Zr);
    element_init(b[i], Zr);
    element_init(m[i], G);
    element_random(a[i]);
    element_random(b[i]);
    element_pow_zn(g0, g, a[i]);
    element_pow_zn(m[i], h, b[i]);
    element_mul(m[i], m[i], g0);
  }

  element_random(asum);
  element_random(bsum);
  element_pow_zn(g0, g, asum);
  element_pow_zn(snark, h, bsum);
  element_mul(snark, snark, g0);

  record(asum, bsum, snark, hole, counter);
  for (;;) {
    int len = element_length_in_bytes(snark);
    unsigned char *buf = pbc_malloc(len);
    unsigned char hash = 0;

    element_to_bytes(buf, snark);
    for (i = 0; i < len; i++) {
      hash += buf[i];
    }
    i = hash % s;
    pbc_free(buf);

    element_mul(snark, snark, m[i]);
    element_add(asum, asum, a[i]);
    element_add(bsum, bsum, b[i]);

    for (i = 0; i < hole->count; i++) {
      snapshot_ptr ss = hole->item[i];
      if (!element_cmp(snark, ss->snark)) {
        element_sub(bsum, bsum, ss->b);
        element_sub(asum, ss->a, asum);
        //answer is x such that x * bsum = asum
        //complications arise if gcd(bsum, r) > 1
        //which can happen if r is not prime
        if (!mpz_probab_prime_p(Zr->order, 10)) {
          mpz_t za, zb, zd, zm;

          mpz_init(za);
          mpz_init(zb);
          mpz_init(zd);
          mpz_init(zm);

          element_to_mpz(za, asum);
          element_to_mpz(zb, bsum);
          mpz_gcd(zd, zb, Zr->order);
          mpz_divexact(zm, Zr->order, zd);
          mpz_divexact(zb, zb, zd);
          //if zd does not divide za there is no solution
          mpz_divexact(za, za, zd);
          mpz_invert(zb, zb, zm);
          mpz_mul(zb, za, zb);
          mpz_mod(zb, zb, zm);
          do {
            element_pow_mpz(g0, g, zb);
            if (!element_cmp(g0, h)) {
              element_set_mpz(x, zb);
              break;
            }
            mpz_add(zb, zb, zm);
            mpz_sub_ui(zd, zd, 1);
          } while (mpz_sgn(zd));
          mpz_clear(zm);
          mpz_clear(za);
          mpz_clear(zb);
          mpz_clear(zd);
        } else {
          element_div(x, asum, bsum);
        }
        found = 1;
        break;
      }
    }
    if (found) break;

    mpz_add_ui(counter, counter, 1);
    if (mpz_tstbit(counter, interval)) {
      record(asum, bsum, snark, hole, counter);
      interval++;
    }
  }

  for (i = 0; i < s; i++) {
    element_clear(a[i]);
    element_clear(b[i]);
    element_clear(m[i]);
  }
  element_clear(g0);
  element_clear(snark);
  for (i = 0; i < hole->count; i++) {
    snapshot_ptr ss = hole->item[i];
    element_clear(ss->a);
    element_clear(ss->b);
    element_clear(ss->snark);
    pbc_free(ss);
  }
  darray_clear(hole);
  element_clear(asum);
  element_clear(bsum);
  mpz_clear(counter);
}
