// Naive implementation of F_p.
// Little more than wrappers around GMP mpz functions.

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

static void zp_init(element_ptr e) {
  e->data = pbc_malloc(sizeof(mpz_t));
  mpz_init(e->data);
}

static void zp_clear(element_ptr e) {
  mpz_clear(e->data);
  pbc_free(e->data);
}

static void zp_set_si(element_ptr e, signed long int op) {
  mpz_set_si(e->data, op);
  mpz_mod(e->data, e->data, e->field->order);
}

static void zp_set_mpz(element_ptr e, mpz_ptr z) {
  mpz_set(e->data, z);
  mpz_mod(e->data, e->data, e->field->order);
}

static void zp_set0(element_ptr e) {
  mpz_set_si(e->data, 0);
}

static void zp_set1(element_ptr e) {
  mpz_set_si(e->data, 1);
}

static size_t zp_out_str(FILE * stream, int base, element_ptr e) {
  return mpz_out_str(stream, base, e->data);
}

static int zp_snprint(char *s, size_t n, element_ptr e) {
  return gmp_snprintf(s, n, "%Zd", e->data);
}

static int zp_set_str(element_ptr e, const char *s, int base) {
  int result = pbc_mpz_set_str(e->data, s, base);
  mpz_mod(e->data, e->data, e->field->order);
  return result;
}

static int zp_sgn_odd(element_ptr a) {
  mpz_ptr z = a->data;

  return mpz_is0(z) ? 0 : (mpz_odd_p(z) ? 1 : -1);
}

static int zp_sgn_even(element_ptr a) {
  mpz_t z;
  mpz_init(z);
  int res;

  if (mpz_is0(a->data)) {
    res = 0;
  } else {
    mpz_add(z, a->data, a->data);
    res = mpz_cmp(z, a->field->order);
  }
  mpz_clear(z);
  return res;
}

static void zp_add(element_ptr n, element_ptr a, element_ptr b) {
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

static void zp_sub(element_ptr n, element_ptr a, element_ptr b) {
  //mpz_sub(n->data, a->data, b->data);
  //mpz_mod(n->data, n->data, n->field->order);
  mpz_sub(n->data, a->data, b->data);
  if (mpz_sgn((mpz_ptr) n->data) < 0) {
    mpz_add(n->data, n->data, n->field->order);
  }
}

static void zp_square(element_ptr c, element_ptr a) {
  /*
     mpz_mul(c->data, a->data, a->data);
     mpz_mod(c->data, c->data, c->field->order);
   */
  mpz_powm_ui(c->data, a->data, 2, c->field->order);

  /*
     const mpz_ptr prime = c->field->order;
     const size_t t = prime->_mp_size;
     const mpz_ptr p = a->data;
     const mpz_ptr r = c->data;
     mp_limb_t tmp[2 * t];
     mp_limb_t qp[t + 1];

     mpn_mul_n(tmp, p->_mp_d, p->_mp_d, t);

     mpn_tdiv_qr(qp, r->_mp_d, 0, tmp, 2 * t, prime->_mp_d, t);
   */
}

static void zp_double(element_ptr n, element_ptr a) {
  //mpz_add(n->data, a->data, a->data);
  mpz_mul_2exp(n->data, a->data, 1);
  if (mpz_cmp(n->data, n->field->order) >= 0) {
    mpz_sub(n->data, n->data, n->field->order);
  }
}

static void zp_halve(element_ptr n, element_ptr a) {
  mpz_ptr z = a->data;
  if (mpz_odd_p(z)) {
    mpz_add(n->data, z, a->field->order);
    mpz_tdiv_q_2exp(n->data, n->data, 1);
  } else {
    mpz_tdiv_q_2exp(n->data, a->data, 1);
  }
}

static void zp_mul(element_ptr n, element_ptr a, element_ptr b) {
  mpz_mul(n->data, a->data, b->data);
  mpz_mod(n->data, n->data, n->field->order);
}

static void zp_mul_mpz(element_ptr n, element_ptr a, mpz_ptr z) {
  mpz_mul(n->data, a->data, z);
  mpz_mod(n->data, n->data, n->field->order);
}

static void zp_mul_si(element_ptr n, element_ptr a, signed long int z) {
  mpz_mul_si(n->data, a->data, z);
  mpz_mod(n->data, n->data, n->field->order);
}

static void zp_pow_mpz(element_ptr n, element_ptr a, mpz_ptr z) {
  mpz_powm(n->data, a->data, z, n->field->order);
}

static void zp_set(element_ptr n, element_ptr a) {
  mpz_set(n->data, a->data);
}

static void zp_neg(element_ptr n, element_ptr a) {
  if (mpz_is0(a->data)) {
    mpz_set_ui(n->data, 0);
  } else {
    mpz_sub(n->data, n->field->order, a->data);
  }
}

static void zp_invert(element_ptr n, element_ptr a) {
  mpz_invert(n->data, a->data, n->field->order);
}

static void zp_random(element_ptr n) {
  pbc_mpz_random(n->data, n->field->order);
}

static void zp_from_hash(element_ptr n, void *data, int len) {
  pbc_mpz_from_hash(n->data, n->field->order, data, len);
}

static int zp_is1(element_ptr n) {
  return !mpz_cmp_ui((mpz_ptr) n->data, 1);
}

static int zp_is0(element_ptr n) {
  return mpz_is0(n->data);
}

static int zp_cmp(element_ptr a, element_ptr b) {
  return mpz_cmp((mpz_ptr) a->data, (mpz_ptr) b->data);
}

static int zp_is_sqr(element_ptr a) {
  //0 is a square
  if (mpz_is0(a->data)) return 1;
  return mpz_legendre(a->data, a->field->order) == 1;
}

static void zp_field_clear(field_t f) {
  UNUSED_VAR(f);
}

static int zp_to_bytes(unsigned char *data, element_t e) {
  int n;

  n = e->field->fixed_length_in_bytes;

  pbc_mpz_out_raw_n(data, n, e->data);
  return n;
}

static int zp_from_bytes(element_t e, unsigned char *data) {
  mpz_ptr z = e->data;
  int n;
  n = e->field->fixed_length_in_bytes;
  mpz_import(z, n, 1, 1, 1, 0, data);
  return n;
}

static void zp_to_mpz(mpz_ptr z, element_ptr a) {
  mpz_set(z, a->data);
}

static void zp_out_info(FILE * out, field_ptr f) {
  element_fprintf(out, "GF(%Zd), GMP wrapped", f->order);
}

void field_init_naive_fp(field_ptr f, mpz_t prime) {
  field_init(f);
  f->init = zp_init;
  f->clear = zp_clear;
  f->set_si = zp_set_si;
  f->set_mpz = zp_set_mpz;
  f->out_str = zp_out_str;
  f->snprint = zp_snprint;
  f->set_str = zp_set_str;
  f->sign = mpz_odd_p(prime) ? zp_sgn_odd : zp_sgn_even;
  f->add = zp_add;
  f->sub = zp_sub;
  f->set = zp_set;
  f->square = zp_square;
  f->doub = zp_double;
  f->halve = zp_halve;
  f->mul = zp_mul;
  f->mul_mpz = zp_mul_mpz;
  f->mul_si = zp_mul_si;
  f->pow_mpz = zp_pow_mpz;
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
  f->sqrt = element_tonelli;
  f->field_clear = zp_field_clear;
  f->to_bytes = zp_to_bytes;
  f->from_bytes = zp_from_bytes;
  f->to_mpz = zp_to_mpz;

  f->out_info = zp_out_info;

  mpz_set(f->order, prime);
  f->data = NULL;
  f->fixed_length_in_bytes = (mpz_sizeinbase(prime, 2) + 7) / 8;
}
