// F_p for small p, i.e. at most sizeof(long) bytes long.
// Assumes long long is at least twice long.

// TODO: Fix outstanding bugs and use in PBC.

#include <stdarg.h>
#include <stdio.h>
#include <stdint.h> // for intptr_t
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_field.h"
#include "pbc_random.h"
#include "pbc_fp.h"
#include "pbc_memory.h"

// Mostly wrappers. We use GMP routines for pow_mpz and invert.

static void fp_init(element_ptr e) {
  unsigned long *p = e->data = pbc_malloc(sizeof(unsigned long));
  *p = 0;
}

static void fp_clear(element_ptr e) {
  pbc_free(e->data);
}

static void fp_set_mpz(element_ptr e, mpz_ptr z) {
  mpz_t r;
  mpz_init(r);
  unsigned long *p = e->field->data;
  unsigned long *l = e->data;
  mpz_fdiv_r_ui(r, z, *p);
  *l = mpz_get_ui(r);
  mpz_clear(r);
}

static void fp_set_si(element_ptr e, signed long int op) {
  unsigned long int *d = e->data;
  unsigned long *p = e->field->data;
  if (op < 0) {
    *d = (-op) % *p;
    *d = *p - *d;
  } else {
    *d = op % *p;
  }
}

static void fp_to_mpz(mpz_ptr z, element_ptr e) {
  unsigned long int *l = e->data;
  mpz_set_ui(z, *l);
}

static void fp_set0(element_ptr e) {
  unsigned long int *l = e->data;
  *l = 0;
}

static void fp_set1(element_ptr e) {
  unsigned long int *l = e->data;
  *l = 1;
}

static int fp_is1(element_ptr e) {
  unsigned long int *l = e->data;
  return *l == 1;
}

static int fp_is0(element_ptr e) {
  unsigned long int *l = e->data;
  return *l == 0;
}

static size_t fp_out_str(FILE *stream, int base, element_ptr e) {
  size_t result;
  mpz_t z;
  mpz_init(z);
  fp_to_mpz(z, e);
  result = mpz_out_str(stream, base, z);
  mpz_clear(z);
  return result;
}

static void fp_add(element_ptr c, element_ptr a, element_ptr b) {
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

static void fp_double(element_ptr c, element_ptr a) {
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

static void fp_sub(element_ptr c, element_ptr a, element_ptr b) {
  unsigned long *prime = a->field->data;
  unsigned long *p = a->data;
  unsigned long *q = b->data;
  unsigned long *r = c->data;

  if (*p >= *q) {
    *r = *p - *q;
  } else {
    *r = *prime - *q + *p;
  }
}

static void fp_mul(element_ptr c, element_ptr a, element_ptr b) {
  unsigned long *prime = a->field->data;
  unsigned long *p = a->data;
  unsigned long *q = b->data;
  unsigned long long ll;
  unsigned long *r = c->data;

  ll = *p * *q;
  *r = ll % *prime;
}

static void fp_square(element_ptr c, element_ptr a) {
  fp_mul(c, a, a);
}

static void fp_neg(element_ptr c, element_ptr a) {
  unsigned long *prime = a->field->data;
  unsigned long *r = c->data;
  unsigned long *p = a->data;
  if (*p) {
    *r = *prime - *p;
  } else {
    *r = 0;
  }
}

static void fp_mul_si(element_ptr c, element_ptr a, signed long int op) {
  unsigned long *prime = a->field->data;
  unsigned long *p = a->data;
  unsigned long long ll;
  unsigned long *r = c->data;

  ll = *p * op;
  *r = ll % *prime;
}

static void fp_pow_mpz(element_ptr c, element_ptr a, mpz_ptr op) {
  unsigned long *r = c->data;
  mpz_t z;
  mpz_init(z);
  fp_to_mpz(z, a);
  mpz_powm(z, z, op, a->field->order);
  *r = mpz_get_ui(z);
  mpz_clear(z);
}

static void fp_set(element_ptr c, element_ptr a) {
  unsigned long *p = a->data;
  unsigned long *r = c->data;
  *r = *p;
}

static void fp_invert(element_ptr c, element_ptr a) {
  unsigned long *r = c->data;
  mpz_t z;
  mpz_init(z);
  fp_to_mpz(z, a);
  mpz_invert(z, z, a->field->order);
  *r = mpz_get_ui(z);
  mpz_clear(z);
}

static void fp_random(element_ptr c) {
  unsigned long *r = c->data;
  mpz_t z;
  mpz_init(z);
  pbc_mpz_random(z, c->field->order);
  *r = mpz_get_ui(z);
  mpz_clear(z);
}

static void fp_from_hash(element_ptr n, void *data, int len) {
  mpz_t z;

  mpz_init(z);
  mpz_import(z, len, -1, 1, -1, 0, data);
  fp_set_mpz(n, z);
  mpz_clear(z);
}

static int fp_cmp(element_ptr a, element_ptr b) {
  unsigned long *p = a->data;
  unsigned long *q = b->data;
  return *p != *q;
}

static int fp_sgn_odd(element_ptr a) {
  unsigned long *p = a->data;
  if (!*p) return 0;
  return *p & 1 ? 1 : -1;
}

static int fp_is_sqr(element_ptr a) {
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

static int fp_to_bytes(unsigned char *data, element_t e) {
  unsigned long *p = e->data;
  unsigned long l = *p;
  int i, n = e->field->fixed_length_in_bytes;
  for (i = 0; i < n; i++) {
    data[n - i - 1] = (unsigned char) l;
    l >>= 8;
  }
  return n;
}

static int fp_from_bytes(element_t e, unsigned char *data) {
  unsigned char *ptr = data;
  unsigned long *p = e->data;
  int i, n = e->field->fixed_length_in_bytes;
  *p = 0;
  for (i=0; i<n; i++) {
    *p <<= 8;
    *p += *ptr;
    ptr++;
  }
  return n;
}

static void fp_field_clear(field_t f) {
  pbc_free(f->data);
}

void field_init_tiny_fp(field_ptr f, mpz_t prime) {
  unsigned long *p;

  PBC_ASSERT(mpz_fits_ulong_p(prime), "modulus too big");

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
  f->sign = fp_sgn_odd;
  f->invert = fp_invert;
  f->random = fp_random;
  f->from_hash = fp_from_hash;
  f->is1 = fp_is1;
  f->is0 = fp_is0;
  f->set0 = fp_set0;
  f->set1 = fp_set1;
  f->is_sqr = fp_is_sqr;
  f->sqrt = element_tonelli;
  f->field_clear = fp_field_clear;
  f->to_bytes = fp_to_bytes;
  f->from_bytes = fp_from_bytes;
  f->to_mpz = fp_to_mpz;

  p = f->data = pbc_malloc(sizeof(long));
  *p = mpz_get_ui(prime);
  {
    unsigned long int l = 255;
    f->fixed_length_in_bytes = 1;
    while (l < *p) {
      f->fixed_length_in_bytes++;
      l <<= 8;
      l += 255;
    }
  }
  mpz_set(f->order, prime);
}
