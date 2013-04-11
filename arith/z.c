// The ring Z.
//
// Wrappers around GMP mpz functions.
#include <stdarg.h>
#include <stdio.h>
#include <stdint.h> // for intptr_t
#include <stdlib.h>
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_field.h"
#include "pbc_z.h"
#include "pbc_random.h"
#include "pbc_fp.h"
#include "pbc_memory.h"

static void z_init(element_ptr e) {
  e->data = pbc_malloc(sizeof(mpz_t));
  mpz_init(e->data);
}

static void z_clear(element_ptr e) {
  mpz_clear(e->data);
  pbc_free(e->data);
}

static void z_set_si(element_ptr e, signed long int op) {
  mpz_set_si(e->data, op);
}

static void z_set_mpz(element_ptr e, mpz_ptr z) {
  mpz_set(e->data, z);
}

static void z_set0(element_ptr e) {
  mpz_set_ui(e->data, 0);
}

static void z_set1(element_ptr e) {
  mpz_set_ui(e->data, 1);
}

static size_t z_out_str(FILE *stream, int base, element_ptr e) {
  return mpz_out_str(stream, base, e->data);
}

static int z_sgn(element_ptr a) {
  mpz_ptr z = a->data;
  return mpz_sgn(z);
}

static void z_add(element_ptr n, element_ptr a, element_ptr b) {
  mpz_add(n->data, a->data, b->data);
}

static void z_sub(element_ptr n, element_ptr a, element_ptr b) {
  mpz_sub(n->data, a->data, b->data);
}

static void z_square(element_ptr c, element_ptr a) {
  mpz_mul(c->data, a->data, a->data);
}

static void z_double(element_ptr n, element_ptr a) {
  mpz_mul_2exp(n->data, a->data, 1);
}

static void z_halve(element_ptr n, element_ptr a) {
  mpz_tdiv_q_2exp(n->data, a->data, -1);
}

static void z_mul(element_ptr n, element_ptr a, element_ptr b) {
  mpz_mul(n->data, a->data, b->data);
}

static void z_mul_mpz(element_ptr n, element_ptr a, mpz_ptr z) {
  mpz_mul(n->data, a->data, z);
}

static void z_mul_si(element_ptr n, element_ptr a, signed long int z) {
  mpz_mul_si(n->data, a->data, z);
}

static void z_pow_mpz(element_ptr n, element_ptr a, mpz_ptr z) {
  mpz_pow_ui(n->data, a->data, mpz_get_ui(z));
}

static void z_set(element_ptr n, element_ptr a) {
  mpz_set(n->data, a->data);
}

static void z_neg(element_ptr n, element_ptr a) {
  mpz_neg(n->data, a->data);
}

static void z_invert(element_ptr n, element_ptr a) {
  if (!mpz_cmpabs_ui(a->data, 1)) {
    mpz_set(n->data, a->data);
  } else mpz_set_ui(n->data, 0);
}

static void z_div(element_ptr c, element_ptr a, element_ptr b) {
  mpz_tdiv_q(c->data, a->data, b->data);
}

//(doesn't make sense if order is infinite)
static void z_random(element_ptr n) {
  mpz_set_ui(n->data, 0);
}

static void z_from_hash(element_ptr n, void *data, int len) {
  mpz_import(n->data, len, -1, 1, -1, 0, data);
}

static int z_is1(element_ptr n) {
  return !mpz_cmp_ui((mpz_ptr) n->data, 1);
}

static int z_is0(element_ptr n) {
  return mpz_is0(n->data);
}

static int z_cmp(element_ptr a, element_ptr b) {
  return mpz_cmp((mpz_ptr) a->data, (mpz_ptr) b->data);
}

static int z_is_sqr(element_ptr a) {
  return mpz_perfect_power_p(a->data);
}

static void z_sqrt(element_ptr c, element_ptr a) {
  mpz_sqrt(c->data, a->data);
}

static void z_field_clear(field_t f) {
  UNUSED_VAR (f);
}

// OpenSSL convention:
//   4 bytes containing length
//   followed by number in big-endian, most-significant bit set if negative
//   (prepending null byte if necessary)
// Positive numbers also the same as mpz_out_raw.
static int z_to_bytes(unsigned char *data, element_t e) {
  mpz_ptr z = e->data;
  size_t msb = mpz_sizeinbase(z, 2);
  size_t n = 4;
  size_t i;

  if (!(msb % 8)) {
    data[4] = 0;
    n++;
  }
  if (mpz_sgn(z) < 0) {
    mpz_export(data + n, NULL, 1, 1, 1, 0, z);
    data[4] |= 128;
  } else {
    mpz_export(data + n, NULL, 1, 1, 1, 0, z);
  }
  n += (msb + 7) / 8 - 4;
  for (i=0; i<4; i++) {
    data[i] = (n >> 8 * (3 - i));
  }
  n += 4;

  return n;
}

static int z_from_bytes(element_t e, unsigned char *data) {
  unsigned char *ptr;
  size_t i, n;
  mpz_ptr z = e->data;
  mpz_t z1;
  int neg = 0;

  mpz_init(z1);
  mpz_set_ui(z, 0);

  ptr = data;
  n = 0;
  for (i=0; i<4; i++) {
    n += ((unsigned int) *ptr) << 8 * (3 - i);
    ptr++;
  }
  if (data[4] & 128) {
    neg = 1;
    data[4] &= 127;
  }
  for (i=0; i<n; i++) {
    mpz_set_ui(z1, *ptr);
    mpz_mul_2exp(z1, z1, 8 * (n - 1 - i));
    ptr++;
    mpz_add(z, z, z1);
  }
  mpz_clear(z1);
  if (neg) mpz_neg(z, z);
  return n;
}

static void z_to_mpz(mpz_ptr z, element_ptr a) {
  mpz_set(z, a->data);
}

static int z_length_in_bytes(element_ptr a) {
  return (mpz_sizeinbase(a->data, 2) + 7) / 8 + 4;
}

static void z_out_info(FILE *out, field_ptr f) {
  UNUSED_VAR(f);
  fprintf(out, "Z: wrapped GMP");
}

static int z_set_str(element_ptr e, const char *s, int base) {
  mpz_t z;
  mpz_init(z);
  int result = pbc_mpz_set_str(z, s, base);
  z_set_mpz(e, z);
  mpz_clear(z);
  return result;
}

void field_init_z(field_ptr f) {
  field_init(f);
  f->init = z_init;
  f->clear = z_clear;
  f->set_si = z_set_si;
  f->set_mpz = z_set_mpz;
  f->set_str = z_set_str;
  f->out_str = z_out_str;
  f->sign = z_sgn;
  f->add = z_add;
  f->sub = z_sub;
  f->set = z_set;
  f->square = z_square;
  f->doub = z_double;
  f->halve = z_halve;
  f->mul = z_mul;
  f->mul_mpz = z_mul_mpz;
  f->mul_si = z_mul_si;
  f->pow_mpz = z_pow_mpz;
  f->neg = z_neg;
  f->cmp = z_cmp;
  f->invert = z_invert;
  f->div = z_div;
  f->random = z_random;
  f->from_hash = z_from_hash;
  f->is1 = z_is1;
  f->is0 = z_is0;
  f->set0 = z_set0;
  f->set1 = z_set1;
  f->is_sqr = z_is_sqr;
  f->sqrt = z_sqrt;
  f->field_clear = z_field_clear;
  f->to_bytes = z_to_bytes;
  f->from_bytes = z_from_bytes;
  f->to_mpz = z_to_mpz;
  f->length_in_bytes = z_length_in_bytes;

  f->out_info = z_out_info;

  mpz_set_ui(f->order, 0);
  f->data = NULL;
  f->fixed_length_in_bytes = -1;
}
