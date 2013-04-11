// F_p using Montgomery representation.
//
// Let b = 256^sizeof(mp_limb_t).
// Let R = b^t be the smallest power of b greater than the modulus p.
// Then x is stored as xR (mod p).
// Addition: same as naive implementation.
// Multipication: Montgomery reduction.
// Code assumes the modulus p is odd.
//
// TODO: mul_2exp(x, p->bytes * 8) could be replaced with
// faster code that messes with GMP internals

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

// Per-field data.
typedef struct {
  size_t limbs;           // Number of limbs per element.
  size_t bytes;           // Number of bytes per element.
  mp_limb_t *primelimbs;  // Points to an array of limbs holding the modulus.
  mp_limb_t negpinv;      // -p^-1 mod b
  mp_limb_t *R;           // R mod p
  mp_limb_t *R3;          // R^3 mod p
} *fptr;

// Per-element data.
typedef struct {
  char flag;     // flag == 0 means the element is zero.
  mp_limb_t *d;  // Otherwise d points to an array holding the element.
} *eptr;

// Copies limbs of z into dst and zeroes any leading limbs, where n is the
// total number of limbs.
// Requires z to have at most n limbs.
static inline void set_limbs(mp_limb_t *dst, mpz_t z, size_t n) {
  size_t count;
  mpz_export(dst, &count, -1, sizeof(mp_limb_t), 0, 0, z);
  memset((void *) (((unsigned char *) dst) + count * sizeof(mp_limb_t)),
         0, (n - count) * sizeof(mp_limb_t));
}

static void fp_init(element_ptr e) {
  fptr p = e->field->data;
  eptr ep = e->data = pbc_malloc(sizeof(*ep));
  ep->flag = 0;
  ep->d = pbc_malloc(p->bytes);
}

static void fp_clear(element_ptr e) {
  eptr ep = e->data;
  pbc_free(ep->d);
  pbc_free(e->data);
}

static void fp_set_mpz(element_ptr e, mpz_ptr z) {
  fptr p = e->field->data;
  eptr ep = e->data;
  if (!mpz_sgn(z)) ep->flag = 0;
  else {
    mpz_t tmp;
    mpz_init(tmp);
    mpz_mul_2exp(tmp, z, p->bytes * 8);
    mpz_mod(tmp, tmp, e->field->order);
    if (!mpz_sgn(tmp)) ep->flag = 0;
    else {
      set_limbs(ep->d, tmp, p->limbs);
      ep->flag = 2;
    }
    mpz_clear(tmp);
  }
}

static void fp_set_si(element_ptr e, signed long int op) {
  fptr p = e->field->data;
  eptr ep = e->data;
  if (!op) ep->flag = 0;
  else {
    mpz_t tmp;
    mpz_init(tmp);
    // TODO: Could be optimized.
    mpz_set_si(tmp, op);
    mpz_mul_2exp(tmp, tmp, p->bytes * 8);
    mpz_mod(tmp, tmp, e->field->order);
    if (!mpz_sgn(tmp)) ep->flag = 0;
    else {
      set_limbs(ep->d, tmp, p->limbs);
      ep->flag = 2;
    }
    mpz_clear(tmp);
  }
}

// Montgomery reduction.
// Algorithm II.4 from Blake, Seroussi and Smart.
static void mont_reduce(mp_limb_t *x, mp_limb_t *y, fptr p) {
  size_t t = p->limbs;
  size_t i;
  mp_limb_t flag = 0;
  for (i = 0; i < t; i++) {
    mp_limb_t u = y[i] * p->negpinv;
    mp_limb_t carry = mpn_addmul_1(&y[i], p->primelimbs, t, u);
    //mpn_add_1(&y[i+t], &y[i+t], t - i + 1, carry);
    flag += mpn_add_1(&y[i + t], &y[i + t], t - i, carry);
  }
  if (flag || mpn_cmp(&y[t], p->primelimbs, t) >= 0) {
    mpn_sub_n(x, &y[t], p->primelimbs, t);
  } else {
    // TODO: GMP set might be faster.
    memcpy(x, &y[t], t * sizeof(mp_limb_t));
  }
}

static void fp_to_mpz(mpz_ptr z, element_ptr e) {
  eptr ep = e->data;
  if (!ep->flag) mpz_set_ui(z, 0);
  else {
    // x is stored as xR.
    // We must divide out R to convert to standard representation.
    fptr p = e->field->data;
    mp_limb_t tmp[2 * p->limbs];

    memcpy(tmp, ep->d, p->limbs * sizeof(mp_limb_t));
    memset(&tmp[p->limbs], 0, p->limbs * sizeof(mp_limb_t));
    _mpz_realloc(z, p->limbs);
    mont_reduce(z->_mp_d, tmp, p);
    // Remove leading zero limbs.
    for (z->_mp_size = p->limbs; !z->_mp_d[z->_mp_size - 1]; z->_mp_size--);
  }
}

static void fp_set0(element_ptr e) {
  eptr ep = e->data;
  ep->flag = 0;
}

static void fp_set1(element_ptr e) {
  fptr p = e->field->data;
  eptr ep = e->data;
  ep->flag = 2;
  memcpy(ep->d, p->R, p->bytes);
}

static int fp_is1(element_ptr e) {
  eptr ep = e->data;
  if (!ep->flag) return 0;
  else {
    fptr p = e->field->data;
    return !mpn_cmp(ep->d, p->R, p->limbs);
  }
}

static int fp_is0(element_ptr e) {
  eptr ep = e->data;
  return !ep->flag;
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

static int fp_snprint(char *s, size_t n, element_ptr e) {
  int result;
  mpz_t z;
  mpz_init(z);
  fp_to_mpz(z, e);
  result = gmp_snprintf(s, n, "%Zd", z);
  mpz_clear(z);
  return result;
}

static int fp_set_str(element_ptr e, const char *s, int base) {
  mpz_t z;
  mpz_init(z);
  int result = pbc_mpz_set_str(z, s, base);
  mpz_mod(z, z, e->field->order);
  fp_set_mpz(e, z);
  mpz_clear(z);
  return result;
}

static void fp_set(element_ptr c, element_ptr a) {
  eptr ad = a->data;
  eptr cd = c->data;
  if (a == c) return;
  if (!ad->flag) cd->flag = 0;
  else {
    fptr p = a->field->data;

    // Assembly is faster, but I don't want to stoop to that level.
    // Instead of memcpy(), we rewrite so GMP assembly ends up being invoked.
    /*
       memcpy(cd->d, ad->d, p->bytes);
     */
    mpz_t z1, z2;
    z1->_mp_d = cd->d;
    z2->_mp_d = ad->d;
    z1->_mp_size = z1->_mp_alloc = z2->_mp_size = z2->_mp_alloc = p->limbs;
    mpz_set(z1, z2);

    cd->flag = 2;
  }
}

static void fp_add(element_ptr c, element_ptr a, element_ptr b) {
  eptr ad = a->data, bd = b->data;

  if (!ad->flag) {
    fp_set(c, b);
  } else if (!bd->flag) {
    fp_set(c, a);
  } else {
    eptr cd = c->data;
    fptr p = a->field->data;
    const size_t t = p->limbs;
    mp_limb_t carry;
    carry = mpn_add_n(cd->d, ad->d, bd->d, t);

    if (carry) {
      // Assumes result of following sub is not zero,
      // i.e. modulus cannot be 2^(n * bits_per_limb).
      mpn_sub_n(cd->d, cd->d, p->primelimbs, t);
      cd->flag = 2;
    } else {
      int i = mpn_cmp(cd->d, p->primelimbs, t);
      if (!i) {
        cd->flag = 0;
      } else {
        cd->flag = 2;
        if (i > 0) {
          mpn_sub_n(cd->d, cd->d, p->primelimbs, t);
        }
      }
    }
  }
}

static void fp_double(element_ptr c, element_ptr a) {
  eptr ad = a->data, cd = c->data;
  if (!ad->flag) {
    cd->flag = 0;
  } else {
    fptr p = c->field->data;
    const size_t t = p->limbs;
    if (mpn_lshift(cd->d, ad->d, t, 1)) {
      cd->flag = 2;
      // Again, assumes result is not zero.
      mpn_sub_n(cd->d, cd->d, p->primelimbs, t);
    } else {
      int i = mpn_cmp(cd->d, p->primelimbs, t);
      if (!i) {
        cd->flag = 0;
      } else {
        cd->flag = 2;
        if (i > 0) {
          mpn_sub_n(cd->d, cd->d, p->primelimbs, t);
        }
      }
    }
  }
}

static void fp_halve(element_ptr c, element_ptr a) {
  eptr ad = a->data, cd = c->data;
  if (!ad->flag) {
    cd->flag = 0;
  } else {
    fptr p = c->field->data;
    const size_t t = p->limbs;
    int carry = 0;
    mp_limb_t *alimb = ad->d;
    mp_limb_t *climb = cd->d;
    if (alimb[0] & 1) {
      carry = mpn_add_n(climb, alimb, p->primelimbs, t);
    } else fp_set(c, a);

    mpn_rshift(climb, climb, t, 1);
    if (carry) climb[t - 1] |= ((mp_limb_t) 1) << (sizeof(mp_limb_t) * 8 - 1);
  }
}

static void fp_neg(element_ptr c, element_ptr a) {
  eptr ad = a->data, cd = c->data;
  if (!ad->flag) cd->flag = 0;
  else {
    fptr p = a->field->data;
    mpn_sub_n(cd->d, p->primelimbs, ad->d, p->limbs);
    cd->flag = 2;
  }
}

static void fp_sub(element_ptr c, element_ptr a, element_ptr b) {
  eptr ad = a->data, bd = b->data;

  if (!ad->flag) {
    fp_neg(c, b);
  } else if (!bd->flag) {
    fp_set(c, a);
  } else {
    fptr p = c->field->data;
    size_t t = p->limbs;
    eptr cd = c->data;
    int i = mpn_cmp(ad->d, bd->d, t);

    if (i == 0) {
      cd->flag = 0;
    } else {
      cd->flag = 2;
      mpn_sub_n(cd->d, ad->d, bd->d, t);
      if (i < 0) {
        mpn_add_n(cd->d, cd->d, p->primelimbs, t);
      }
    }
  }
}

// Montgomery multiplication.
// See Blake, Seroussi and Smart.
static inline void mont_mul(mp_limb_t *c, mp_limb_t *a, mp_limb_t *b,
                            fptr p) {
  // Instead of right shifting every iteration
  // I allocate more room for the z array.
  size_t i, t = p->limbs;
  mp_limb_t z[2 * t + 1];
  mp_limb_t u = (a[0] * b[0]) * p->negpinv;
  mp_limb_t v = z[t] = mpn_mul_1(z, b, t, a[0]);
  z[t] += mpn_addmul_1(z, p->primelimbs, t, u);
  z[t + 1] = z[t] < v;  // Handle overflow.
  for (i = 1; i < t; i++) {
    u = (z[i] + a[i] * b[0]) * p->negpinv;
    v = z[t + i] += mpn_addmul_1(z + i, b, t, a[i]);
    z[t + i] += mpn_addmul_1(z + i, p->primelimbs, t, u);
    z[t + i + 1] = z[t + i] < v;
  }
  if (z[t * 2] || mpn_cmp(z + t, p->primelimbs, t) >= 0) {
    mpn_sub_n(c, z + t, p->primelimbs, t);
  } else {
    memcpy(c, z + t, t * sizeof(mp_limb_t));
    // Doesn't seem to make a difference:
    /*
       mpz_t z1, z2;
       z1->_mp_d = c;
       z2->_mp_d = z + t;
       z1->_mp_size = z1->_mp_alloc = z2->_mp_size = z2->_mp_alloc = t;
       mpz_set(z1, z2);
     */
  }
}

static void fp_mul(element_ptr c, element_ptr a, element_ptr b) {
  eptr ad = a->data, bd = b->data;
  eptr cd = c->data;

  if (!ad->flag || !bd->flag) {
    cd->flag = 0;
  } else {
    fptr p = c->field->data;
    mont_mul(cd->d, ad->d, bd->d, p);
    cd->flag = 2;
  }
}

static void fp_pow_mpz(element_ptr c, element_ptr a, mpz_ptr op) {
  // Alternative: rewrite GMP mpz_powm().
  fptr p = a->field->data;
  eptr ad = a->data;
  eptr cd = c->data;
  if (!ad->flag) cd->flag = 0;
  else {
    mpz_t z;
    mpz_init(z);
    fp_to_mpz(z, a);
    mpz_powm(z, z, op, a->field->order);
    mpz_mul_2exp(z, z, p->bytes * 8);
    mpz_mod(z, z, a->field->order);
    set_limbs(cd->d, z, p->limbs);
    mpz_clear(z);
    cd->flag = 2;
  }
}

// Inversion is slower than in a naive Fp implementation because of an extra
// multiplication.
// Requires nonzero a.
static void fp_invert(element_ptr c, element_ptr a) {
  eptr ad = a->data;
  eptr cd = c->data;
  fptr p = a->field->data;
  mp_limb_t tmp[p->limbs];
  mpz_t z;

  mpz_init(z);

  // Copy the limbs into a regular mpz_t so we can invert using the standard
  // mpz_invert().
  mpz_import(z, p->limbs, -1, sizeof(mp_limb_t), 0, 0, ad->d);
  mpz_invert(z, z, a->field->order);
  set_limbs(tmp, z, p->limbs);

  // Normalize.
  mont_mul(cd->d, tmp, p->R3, p);
  cd->flag = 2;
  mpz_clear(z);
}

static void fp_random(element_ptr a) {
  fptr p = a->field->data;
  eptr ad = a->data;
  mpz_t z;
  mpz_init(z);
  pbc_mpz_random(z, a->field->order);
  if (mpz_sgn(z)) {
    mpz_mul_2exp(z, z, p->bytes * 8);
    mpz_mod(z, z, a->field->order);
    set_limbs(ad->d, z, p->limbs);
    ad->flag = 2;
  } else {
    ad->flag = 0;
  }
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
  eptr ad = a->data, bd = b->data;
  if (!ad->flag) return bd->flag;
  else {
    fptr p = a->field->data;
    return mpn_cmp(ad->d, bd->d, p->limbs);
    //return memcmp(ad->d, bd->d, p->limbs);
  }
}

static int fp_sgn_odd(element_ptr a) {
  eptr ad = a->data;
  if (!ad->flag) return 0;
  else {
    mpz_t z;
    mpz_init(z);
    int res;
    fp_to_mpz(z, a);
    res = mpz_odd_p(z) ? 1 : -1;
    mpz_clear(z);
    return res;
  }
}

static int fp_is_sqr(element_ptr a) {
  eptr ad = a->data;
  int res;
  mpz_t z;
  mpz_init(z);
  // 0 is a square.
  if (!ad->flag) return 1;
  fp_to_mpz(z, a);
  res = mpz_legendre(z, a->field->order) == 1;
  mpz_clear(z);
  return res;
}

static int fp_to_bytes(unsigned char *data, element_t a) {
  mpz_t z;
  int n = a->field->fixed_length_in_bytes;

  mpz_init(z);
  fp_to_mpz(z, a);
  pbc_mpz_out_raw_n(data, n, z);
  mpz_clear(z);
  return n;
}

static int fp_from_bytes(element_t a, unsigned char *data) {
  fptr p = a->field->data;
  eptr ad = a->data;
  int n;
  mpz_t z;

  mpz_init(z);

  n = a->field->fixed_length_in_bytes;
  mpz_import(z, n, 1, 1, 1, 0, data);
  if (!mpz_sgn(z)) ad->flag = 0;
  else {
    ad->flag = 2;
    mpz_mul_2exp(z, z, p->bytes * 8);
    mpz_mod(z, z, a->field->order);
    set_limbs(ad->d, z, p->limbs);
  }
  mpz_clear(z);
  return n;
}

static void fp_field_clear(field_t f) {
  fptr p = f->data;
  pbc_free(p->primelimbs);
  pbc_free(p->R);
  pbc_free(p->R3);
  pbc_free(p);
}

// The only public functions. All the above should be static.

static void fp_out_info(FILE * out, field_ptr f) {
  element_fprintf(out, "GF(%Zd): Montgomery representation", f->order);
}

void field_init_mont_fp(field_ptr f, mpz_t prime) {
  PBC_ASSERT(!mpz_fits_ulong_p(prime), "modulus too small");
  fptr p;
  field_init(f);
  f->init = fp_init;
  f->clear = fp_clear;
  f->set_si = fp_set_si;
  f->set_mpz = fp_set_mpz;
  f->out_str = fp_out_str;
  f->snprint = fp_snprint;
  f->set_str = fp_set_str;
  f->add = fp_add;
  f->sub = fp_sub;
  f->set = fp_set;
  f->mul = fp_mul;
  f->doub = fp_double;
  f->halve = fp_halve;
  f->pow_mpz = fp_pow_mpz;
  f->neg = fp_neg;
  f->sign = fp_sgn_odd;
  f->cmp = fp_cmp;
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
  f->out_info = fp_out_info;

  // Initialize per-field data specific to this implementation.
  p = f->data = pbc_malloc(sizeof(*p));
  p->limbs = mpz_size(prime);
  p->bytes = p->limbs * sizeof(mp_limb_t);
  p->primelimbs = pbc_malloc(p->bytes);
  mpz_export(p->primelimbs, &p->limbs, -1, sizeof(mp_limb_t), 0, 0, prime);

  mpz_set(f->order, prime);
  f->fixed_length_in_bytes = (mpz_sizeinbase(prime, 2) + 7) / 8;

  // Compute R, R3 and negpinv.
  mpz_t z;
  mpz_init(z);

  p->R = pbc_malloc(p->bytes);
  p->R3 = pbc_malloc(p->bytes);
  mpz_setbit(z, p->bytes * 8);
  mpz_mod(z, z, prime);
  set_limbs(p->R, z, p->limbs);

  mpz_powm_ui(z, z, 3, prime);
  set_limbs(p->R3, z, p->limbs);

  mpz_set_ui(z, 0);

  // Algorithm II.5 in Blake, Seroussi and Smart is better but this suffices
  // since we're only doing it once.
  mpz_setbit(z, p->bytes * 8);
  mpz_invert(z, prime, z);
  p->negpinv = -mpz_get_ui(z);
  mpz_clear(z);
}
