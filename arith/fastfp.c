// Naive implementation of F_p.
// Uses lowlevel GMP routines (mpn_* functions).
//
// Within an element_t, ''data'' field of element holds pointer to array of
// mp_limb_t, which is allocated on init and freed on clear.
// Its size is fixed and determined by the number of limbs in the modulus.
// This simplifies code but is inefficient for storing values like 0 and 1.

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

struct fp_field_data_s {
  size_t limbs;
  size_t bytes;
  mp_limb_t *primelimbs;
};
typedef struct fp_field_data_s fp_field_data_t[1];
typedef struct fp_field_data_s *fp_field_data_ptr;

static void fp_init(element_ptr e) {
  fp_field_data_ptr p = e->field->data;
  e->data = pbc_malloc(p->bytes);
  memset(e->data, 0, p->bytes);
  // e->data = pbc_calloc(sizeof(mp_limb_t), p->limbs);
}

static void fp_clear(element_ptr e) {
  pbc_free(e->data);
}

static inline void from_mpz(element_ptr e, mpz_ptr z) {
  fp_field_data_ptr p = e->field->data;
  size_t count;
  mpz_export(e->data, &count, -1, sizeof(mp_limb_t), 0, 0, z);
  memset((void *) (((unsigned char *) e->data) + count * sizeof(mp_limb_t)), 0,
         (p->limbs - count) * sizeof(mp_limb_t));
}

static void fp_set_mpz(element_ptr e, mpz_ptr z) {
  mpz_t tmp;
  mpz_init(tmp);
  mpz_mod(tmp, z, e->field->order);
  from_mpz(e, tmp);
  mpz_clear(tmp);
}

static void fp_set_si(element_ptr e, signed long int op) {
  const fp_field_data_ptr p = e->field->data;
  const size_t t = p->limbs;
  mp_limb_t *d = e->data;
  if (op < 0) {
    mpn_sub_1(d, p->primelimbs, t, -op);
  } else {
    d[0] = op;
    memset(&d[1], 0, sizeof(mp_limb_t) * (t - 1));
  }
}

static void fp_to_mpz(mpz_ptr z, element_ptr a) {
  fp_field_data_ptr p = a->field->data;
  mpz_import(z, p->limbs, -1, sizeof(mp_limb_t), 0, 0, a->data);
}

static void fp_set0(element_ptr e) {
  fp_field_data_ptr p = e->field->data;
  memset(e->data, 0, p->bytes);
}

static void fp_set1(element_ptr e) {
  fp_field_data_ptr p = e->field->data;
  mp_limb_t *d = e->data;
  memset(&d[1], 0, p->bytes - sizeof(mp_limb_t));
  d[0] = 1;
}

static int fp_is1(element_ptr e) {
  fp_field_data_ptr p = e->field->data;
  size_t i, t = p->limbs;
  mp_limb_t *d = e->data;
  if (d[0] != 1) return 0;
  for (i = 1; i < t; i++) if (d[i]) return 0;
  return 1;
}

static int fp_is0(element_ptr e) {
  fp_field_data_ptr p = e->field->data;
  size_t i, t = p->limbs;
  mp_limb_t *d = e->data;
  for (i = 0; i < t; i++) if (d[i]) return 0;
  return 1;
}

static size_t fp_out_str(FILE * stream, int base, element_ptr e) {
  size_t result;
  mpz_t z;
  mpz_init(z);
  fp_to_mpz(z, e);
  result = mpz_out_str(stream, base, z);
  mpz_clear(z);
  return result;
}

static void fp_add(element_ptr r, element_ptr a, element_ptr b) {
  fp_field_data_ptr p = r->field->data;
  const size_t t = p->limbs;
  mp_limb_t carry;
  carry = mpn_add_n(r->data, a->data, b->data, t);

  if (carry || mpn_cmp(r->data, p->primelimbs, t) >= 0) {
    mpn_sub_n(r->data, r->data, p->primelimbs, t);
  }
}

static void fp_double(element_ptr r, element_ptr a) {
  fp_field_data_ptr p = r->field->data;
  const size_t t = p->limbs;
  if (mpn_lshift(r->data, a->data, t, 1)
      || mpn_cmp(r->data, p->primelimbs, t) >= 0) {
    mpn_sub_n(r->data, r->data, p->primelimbs, t);
  }
}

static void fp_set(element_ptr c, element_ptr a) {
  fp_field_data_ptr p = a->field->data;
  if (c == a) return;

  // Assembly is faster here, but I don't want to stoop to that level.
  // Instead of calling slower memcpy, wrap stuff so that GMP assembly
  // gets called.
  /*
     memcpy(c->data, a->data, p->bytes);
   */
  mpz_t z1, z2;
  z1->_mp_d = c->data;
  z2->_mp_d = a->data;
  z1->_mp_size = z1->_mp_alloc = z2->_mp_size = z2->_mp_alloc = p->limbs;
  mpz_set(z1, z2);
}

static void fp_halve(element_ptr r, element_ptr a) {
  fp_field_data_ptr p = r->field->data;
  const size_t t = p->limbs;
  int carry = 0;
  mp_limb_t *alimb = a->data;
  mp_limb_t *rlimb = r->data;
  if (alimb[0] & 1) carry = mpn_add_n(rlimb, alimb, p->primelimbs, t);
  else fp_set(r, a);

  mpn_rshift(rlimb, rlimb, t, 1);
  if (carry) rlimb[t - 1] |= ((mp_limb_t) 1) << (sizeof(mp_limb_t) * 8 - 1);
}

static void fp_sub(element_ptr r, element_ptr a, element_ptr b) {
  fp_field_data_ptr p = r->field->data;
  size_t t = p->limbs;
  if (mpn_sub_n(r->data, a->data, b->data, t)) {
    mpn_add_n(r->data, r->data, p->primelimbs, t);
  }
}

static void fp_mul(element_ptr c, element_ptr a, element_ptr b) {
  fp_field_data_ptr p = c->field->data;
  size_t t = p->limbs;
  //mp_limb_t tmp[3 * t + 1];
  //mp_limb_t *qp = &tmp[2 * t];
  mp_limb_t tmp[2 * t];
  mp_limb_t qp[t + 1];
  //static mp_limb_t tmp[2 * 100];
  //static mp_limb_t qp[100 + 1];

  mpn_mul_n(tmp, a->data, b->data, t);

  mpn_tdiv_qr(qp, c->data, 0, tmp, 2 * t, p->primelimbs, t);
}

static void fp_square(element_ptr c, element_ptr a) {
  const fp_field_data_ptr r = c->field->data;
  mpz_t z1, z2;
  size_t diff;

  z1->_mp_d = c->data;
  z1->_mp_size = z1->_mp_alloc = r->limbs;
  if (c == a) {
    mpz_powm_ui(z1, z1, 2, c->field->order);
  } else {
    z2->_mp_d = a->data;
    z2->_mp_size = z2->_mp_alloc = r->limbs;
    mpz_powm_ui(z1, z2, 2, c->field->order);
  }

  diff = r->limbs - z1->_mp_size;
  if (diff) memset(&z1->_mp_d[z1->_mp_size], 0, diff * sizeof(mp_limb_t));

  //mpn_sqr_n() might make the code below faster than the code above
  //but GMP doesn't expose this function
  /*
     const fp_field_data_ptr r = c->field->data;
     const size_t t = r->limbs;
     mp_limb_t tmp[2 * t];
     mp_limb_t qp[t + 1];

     mpn_mul_n(tmp, a->data, a->data, t);

     mpn_tdiv_qr(qp, c->data, 0, tmp, 2 * t, r->primelimbs, t);
   */
}

static void fp_neg(element_ptr n, element_ptr a) {
  if (fp_is0(a)) {
    fp_set0(n);
  } else {
    fp_field_data_ptr p = a->field->data;
    mpn_sub_n(n->data, p->primelimbs, a->data, p->limbs);
  }
}

static void fp_mul_si(element_ptr e, element_ptr a, signed long int op) {
  fp_field_data_ptr p = e->field->data;
  size_t t = p->limbs;
  mp_limb_t tmp[t + 1];
  mp_limb_t qp[2];

  tmp[t] = mpn_mul_1(tmp, a->data, t, labs(op));
  mpn_tdiv_qr(qp, e->data, 0, tmp, t + 1, p->primelimbs, t);
  if (op < 0) {
    fp_neg(e, e);
  }
}

static void fp_pow_mpz(element_ptr c, element_ptr a, mpz_ptr op) {
  mpz_t z;
  mpz_init(z);
  fp_to_mpz(z, a);
  mpz_powm(z, z, op, c->field->order);
  from_mpz(c, z);
  mpz_clear(z);
}

static void fp_invert(element_ptr e, element_ptr a) {
  mpz_t z;
  mpz_init(z);
  fp_to_mpz(z, a);
  mpz_invert(z, z, e->field->order);
  from_mpz(e, z);
  mpz_clear(z);
}

static void fp_random(element_ptr a) {
  mpz_t z;
  mpz_init(z);
  pbc_mpz_random(z, a->field->order);
  from_mpz(a, z);
  mpz_clear(z);
}

static void fp_from_hash(element_ptr a, void *data, int len) {
  mpz_t z;

  mpz_init(z);
  pbc_mpz_from_hash(z, a->field->order, data, len);
  fp_set_mpz(a, z);
  mpz_clear(z);
}

static int fp_cmp(element_ptr a, element_ptr b) {
  fp_field_data_ptr p = a->field->data;
  return mpn_cmp(a->data, b->data, p->limbs);
  //return memcmp(a->data, b->data, p->limbs);
}

static int fp_sgn_odd(element_ptr a) {
  if (fp_is0(a)) return 0;
  mp_limb_t *lp = a->data;
  return lp[0] & 1 ? 1 : -1;
}

static int fp_sgn_even(element_ptr a) {
  fp_field_data_ptr p = a->field->data;
  if (fp_is0(a)) return 0;
  mp_limb_t sum[p->limbs];

  int carry = mpn_add_n(sum, a->data, a->data, p->limbs);
  if (carry) return 1;
  return mpn_cmp(sum, p->primelimbs, p->limbs);
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
  mpz_t z;
  int n;

  mpz_init(z);
  fp_to_mpz(z, e);
  n = e->field->fixed_length_in_bytes;
  pbc_mpz_out_raw_n(data, n, z);
  mpz_clear(z);
  return n;
}

static int fp_from_bytes(element_t e, unsigned char *data) {
  int n;
  mpz_t z;

  mpz_init(z);

  n = e->field->fixed_length_in_bytes;
  mpz_import(z, n, 1, 1, 1, 0, data);
  fp_set_mpz(e, z);
  mpz_clear(z);
  return n;
}

static void fp_field_clear(field_t f) {
  fp_field_data_ptr p = f->data;
  pbc_free(p->primelimbs);
  pbc_free(p);
}

void field_init_fast_fp(field_ptr f, mpz_t prime) {
  PBC_ASSERT(!mpz_fits_ulong_p(prime), "modulus too small");
  fp_field_data_ptr p;
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
  f->halve = fp_halve;
  f->pow_mpz = fp_pow_mpz;
  f->neg = fp_neg;
  f->cmp = fp_cmp;
  f->sign = mpz_odd_p(prime) ? fp_sgn_odd : fp_sgn_even;
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

  p = f->data = pbc_malloc(sizeof(fp_field_data_t));
  p->limbs = mpz_size(prime);
  p->bytes = p->limbs * sizeof(mp_limb_t);
  p->primelimbs = pbc_malloc(p->bytes);
  mpz_export(p->primelimbs, &p->limbs, -1, sizeof(mp_limb_t), 0, 0, prime);

  mpz_set(f->order, prime);
  f->fixed_length_in_bytes = (mpz_sizeinbase(prime, 2) + 7) / 8;
}
