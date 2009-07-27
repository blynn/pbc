// Routines common to all implementations of F_p.
//
// The other routines are in specific implementations of F_p: naivefp.c,
// fastfp.c, fasterfp.c and montfp.c. For pairing-based cryptosystems, montfp.c
// is the fastest, but I keep all versions around for testing, and also to
// show off the modularity of the code.

#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <gmp.h>
#include <string.h>
#include "pbc_utils.h"
#include "pbc_field.h"
#include "pbc_fp.h"
#include "pbc_memory.h"

// By default, use the montfp.c implementation of F_p. After
// pbc_tweak_use_fp(), future field_init_fp calls will use the specified
// implementation. This is useful for benchmarking and testing.
static void (*option_fpinit) (field_ptr f, mpz_t prime) = field_init_mont_fp;

void pbc_tweak_use_fp(char *s) {
  if (!strcmp(s, "naive")) {
    option_fpinit = field_init_naive_fp;
  } else if (!strcmp(s, "fast")) {
    option_fpinit = field_init_fast_fp;
  } else if (!strcmp(s, "faster")) {
    option_fpinit = field_init_faster_fp;
  } else if (!strcmp(s, "mont")) {
    option_fpinit = field_init_mont_fp;
  } else {
    pbc_error("no such Fp implementation: %s", s);
  }
}

void field_init_fp(field_ptr f, mpz_t modulus) {
  if (mpz_fits_ulong_p(modulus)) {
    // If this case mattered, I'd have written a F_p implementation specialized
    // for moduli that fits into machine words.
    field_init_naive_fp(f, modulus);
  } else {
    if (mpz_odd_p(modulus)) {
      option_fpinit(f, modulus);
    } else {
      // montfp.c only supports odd moduli.
      field_init_faster_fp(f, modulus);
    }
  }
}

// Square root algorithm.
void fp_tonelli(element_ptr x, element_ptr a) {
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
  for (i = 2; i <= s; i++) {
    mpz_sub_ui(t0, a->field->order, 1);
    mpz_tdiv_q_2exp(t0, t0, i);
    element_pow_mpz(e0, ginv, e);
    element_mul(e0, e0, a);
    element_pow_mpz(e0, e0, t0);
    if (!element_is1(e0)) mpz_setbit(e, i - 1);
  }
  element_pow_mpz(e0, ginv, e);
  element_mul(e0, e0, a);
  mpz_add_ui(t, t, 1);
  mpz_tdiv_q_2exp(t, t, 1);
  mpz_tdiv_q_2exp(e, e, 1);

  // (suggested by Hovav Shacham) replace next three lines with
  //   element_pow2_mpz(x, e0, t, nqr, e);
  // once sliding windows are implemented for pow2
  element_pow_mpz(e0, e0, t);
  element_pow_mpz(x, nqr, e);
  element_mul(x, x, e0);

  mpz_clear(t);
  mpz_clear(e);
  mpz_clear(t0);
  element_clear(ginv);
  element_clear(e0);
}

// Like mpz_set_str except returns number of bytes read and allows trailing
// junk. This simplifies code for parsing elements like "[123, 456]".
// TODO: Handle 0x, 0X and 0 conventions for hexadecimal and octal.
int pbc_mpz_set_str(mpz_t z, char *s, int base) {
  int b, i = 0;
  mpz_set_ui(z, 0);
  if (!base) b = 10;
  else if (base < 2 || base > 36) return 0;
  else b = base;

  for (;;) {
    int j;
    char c = s[i];
    if (!c) break;
    if (isspace(c)) {
      i++;
      continue;
    }
    if (isdigit(c)) {
      j = c - '0';
    } else if (c >= 'A' && c <= 'Z') {
      j = c - 'A';
    } else if (c >= 'a' && c <= 'z') {
      j = c - 'a';
    } else break;

    if (j >= b) break;

    mpz_mul_ui(z, z, b);
    mpz_add_ui(z, z, j);
    i++;
  }
  return i;
}

int pbc_trial_divide(int (*fun)(mpz_t factor, unsigned int multiplicity),
    mpz_t n, mpz_ptr limit) {
  mpz_t p, m;
  mpz_ptr fac;
  unsigned int mul;

  mpz_init(p);
  mpz_init(m);
  mpz_set(m ,n);
  mpz_set_ui(p, 2);

  while (mpz_cmp_ui(m, 1)) {
    if (mpz_probab_prime_p(m, 10)) {
      mpz_set(p, m);
    }
    if (limit && mpz_cmp(p, limit) > 0) {
      mpz_set(p, m);
    }
    if (mpz_divisible_p(m, p)) {
      fac = pbc_malloc(sizeof(mpz_t));
      mul = 0;
      mpz_init(fac);
      mpz_set(fac, p);
      do {
        mpz_divexact(m, m, p);
        mul++;
      } while (mpz_divisible_p(m, p));
      if (fun(fac, mul)) {
        mpz_clear(m);
        mpz_clear(p);
        return 1;
      }
    }
    mpz_nextprime(p, p);
  }

  mpz_clear(m);
  mpz_clear(p);
  return 0;
}

// For each digit of 'n', call fun(). If it returns 1, then return 1 and
// abort. Otherwise return 0.
int pbc_mpz_trickle(int (*fun)(char), int base, mpz_t n) {
  // TODO: Support different bases.
  if (!base) base = 10;
  if (base < 2 || base > 10) {
    pbc_warn("only bases 2 to 10 supported");
    return 1;
  }
  mpz_t d, z, q;
  mpz_init(d);
  mpz_init(z);
  mpz_init(q);
  mpz_set(z, n);
  int res;
  int len;
  mpz_ui_pow_ui(d, base, len = mpz_sizeinbase(z, base));
  if (mpz_cmp(d, z) > 0) {
    len--;
    mpz_divexact_ui(d, d, base);
  }
  while (mpz_cmp_ui(z, base) >= 0) {
    mpz_fdiv_qr(q, z, z, d);
    res = fun('0' + mpz_get_ui(q));
    if (res) goto clean;
    mpz_divexact_ui(d, d, base);
    len--;
  }
  while (len) {
    res = fun('0');
    if (res) goto clean;
    len--;
  }
  res = fun('0' + mpz_get_ui(z));
clean:
  mpz_clear(q);
  mpz_clear(z);
  mpz_clear(d);
  return res;
}
