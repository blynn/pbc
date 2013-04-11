// Naive implementation of F_p.
// It uses lowlevel GMP routines (mpn_* functions) like fastfp.c, but also
// has a flag for the value 0, avoiding many memsets.
//
// I'm thinking of using the flag to also represent 1, -1,
// but that complicates the logic even more, and I believe I need more
// control than GMP is willing to give in order to avoid expensive
// checks for 1, -1 everywhere.
//
// NOTE: does not work for moduli of the form:
//   2^(something * 8 * sizeof(mp_limb_t))
// See comments in add, double code.
// (This kind of integer mod ring deserves its own implementation anyway.)

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

struct data_s {
  int flag;
  mp_limb_t *d;
};
typedef struct data_s *dataptr;

static void fp_init(element_ptr e) {
  fp_field_data_ptr p = e->field->data;
  dataptr dp = e->data = pbc_malloc(sizeof(struct data_s));
  dp->flag = 0;
  dp->d = pbc_malloc(p->bytes);
}

static void fp_clear(element_ptr e) {
  dataptr dp = e->data;
  pbc_free(dp->d);
  pbc_free(e->data);
}

//assumes z is nonzero
static inline void from_mpz(element_ptr e, mpz_ptr z) {
  fp_field_data_ptr p = e->field->data;
  size_t count;
  dataptr dp = e->data;
  mpz_export(dp->d, &count, -1, sizeof(mp_limb_t), 0, 0, z);
  memset((void *) (((unsigned char *) dp->d) + count * sizeof(mp_limb_t)),
         0, (p->limbs - count) * sizeof(mp_limb_t));
}

static void fp_set_mpz(element_ptr e, mpz_ptr z) {
  dataptr dp = e->data;
  if (!mpz_sgn(z)) {
    dp->flag = 0;
  } else {
    mpz_t tmp;
    mpz_init(tmp);
    mpz_mod(tmp, z, e->field->order);
    from_mpz(e, tmp);
    mpz_clear(tmp);
    dp->flag = 2;
  }
}

static void fp_set_si(element_ptr e, signed long int op) {
  dataptr dp = e->data;
  if (!op) {
    dp->flag = 0;
  } else {
    const fp_field_data_ptr p = e->field->data;
    const size_t t = p->limbs;
    if (op < 0) {
      mpn_sub_1(dp->d, p->primelimbs, t, -op);
    } else {
      dp->d[0] = op;
      memset(&dp->d[1], 0, sizeof(mp_limb_t) * (t - 1));
    }
    dp->flag = 2;
  }
}

static void fp_to_mpz(mpz_ptr z, element_ptr e) {
  dataptr dp = e->data;
  if (!dp->flag) {
    mpz_set_ui(z, 0);
  } else {
    fp_field_data_ptr p = e->field->data;
    mpz_import(z, p->limbs, -1, sizeof(mp_limb_t), 0, 0, dp->d);
  }
}

static void fp_set0(element_ptr e) {
  dataptr dp = e->data;
  dp->flag = 0;
}

static void fp_set1(element_ptr e) {
  fp_field_data_ptr p = e->field->data;
  dataptr dp = e->data;
  dp->flag = 2;
  memset(&dp->d[1], 0, p->bytes - sizeof(mp_limb_t));
  dp->d[0] = 1;
}

static int fp_is1(element_ptr e) {
  dataptr dp = e->data;
  if (!dp->flag) return 0;
  else {
    fp_field_data_ptr p = e->field->data;
    size_t i, t = p->limbs;
    if (dp->d[0] != 1) return 0;
    for (i = 1; i < t; i++) if (dp->d[i]) return 0;
    return 1;
  }
}

static int fp_is0(element_ptr e) {
  dataptr dp = e->data;
  return !dp->flag;
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

static void fp_set(element_ptr c, element_ptr a) {
  dataptr ad = a->data;
  dataptr cd = c->data;
  if (a == c) return;
  if (!ad->flag) {
    cd->flag = 0;
  } else {
    fp_field_data_ptr p = a->field->data;

    //Assembly is faster here, but I don't want to stoop to that level.
    //Instead of calling slower memcpy, wrap stuff so that GMP assembly
    //gets called.
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
  dataptr ad = a->data, bd = b->data;

  if (!ad->flag) {
    fp_set(c, b);
  } else if (!bd->flag) {
    fp_set(c, a);
  } else {
    dataptr cd = c->data;
    fp_field_data_ptr p = a->field->data;
    const size_t t = p->limbs;
    mp_limb_t carry;
    carry = mpn_add_n(cd->d, ad->d, bd->d, t);

    if (carry) {
      //assumes result of following sub is not zero,
      //i.e. modulus cannot be 2^(n * bits_per_limb)
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
  dataptr ad = a->data, cd = c->data;
  if (!ad->flag) {
    cd->flag = 0;
  } else {
    fp_field_data_ptr p = c->field->data;
    const size_t t = p->limbs;
    if (mpn_lshift(cd->d, ad->d, t, 1)) {
      cd->flag = 2;
      //again, assumes result is not zero:
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
  dataptr ad = a->data, cd = c->data;
  if (!ad->flag) {
    cd->flag = 0;
  } else {
    fp_field_data_ptr p = c->field->data;
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
  dataptr ad = a->data, cd = c->data;
  if (!ad->flag) cd->flag = 0;
  else {
    fp_field_data_ptr p = a->field->data;
    mpn_sub_n(cd->d, p->primelimbs, ad->d, p->limbs);
    cd->flag = 2;
  }
}

static void fp_sub(element_ptr c, element_ptr a, element_ptr b) {
  dataptr ad = a->data, bd = b->data;

  if (!ad->flag) {
    fp_neg(c, b);
  } else if (!bd->flag) {
    fp_set(c, a);
  } else {
    fp_field_data_ptr p = c->field->data;
    size_t t = p->limbs;
    dataptr cd = c->data;
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

static void fp_mul(element_ptr c, element_ptr a, element_ptr b) {
  dataptr ad = a->data, bd = b->data;
  dataptr cd = c->data;

  if (!ad->flag || !bd->flag) {
    cd->flag = 0;
  } else {
    fp_field_data_ptr p = c->field->data;
    size_t t = p->limbs;
    //mp_limb_t tmp[3 * t + 1];
    //mp_limb_t *qp = &tmp[2 * t];
    mp_limb_t tmp[2 * t];
    mp_limb_t qp[t + 1];
    //static mp_limb_t tmp[2 * 100];
    //static mp_limb_t qp[100 + 1];

    mpn_mul_n(tmp, ad->d, bd->d, t);

    mpn_tdiv_qr(qp, cd->d, 0, tmp, 2 * t, p->primelimbs, t);
    cd->flag = 2;
  }
}

static void fp_square(element_ptr c, element_ptr a) {
  const fp_field_data_ptr p = c->field->data;
  mpz_t z1, z2;
  size_t diff;
  dataptr ad = a->data;
  dataptr cd = c->data;

  if (!ad->flag) {
    cd->flag = 0;
  } else {
    cd->flag = 2;
    z1->_mp_d = cd->d;
    z1->_mp_size = z1->_mp_alloc = p->limbs;
    if (c == a) {
      mpz_powm_ui(z1, z1, 2, c->field->order);
    } else {
      z2->_mp_d = ad->d;
      z2->_mp_size = z2->_mp_alloc = p->limbs;
      mpz_powm_ui(z1, z2, 2, c->field->order);
    }

    diff = p->limbs - z1->_mp_size;
    if (diff) memset(&z1->_mp_d[z1->_mp_size], 0, diff * sizeof(mp_limb_t));

    //mpn_sqr_n() might make the code below faster than the code above
    //but GMP doesn't expose this function
    /*
       const fp_field_data_ptr p = c->field->data;
       const size_t t = p->limbs;
       mp_limb_t tmp[2 * t];
       mp_limb_t qp[t + 1];

       mpn_mul_n(tmp, ad->d, ad->d, t);

       mpn_tdiv_qr(qp, cd->d, 0, tmp, 2 * t, p->primelimbs, t);
     */
  }
}

static void fp_mul_si(element_ptr c, element_ptr a, signed long int op) {
  dataptr ad = a->data;
  dataptr cd = c->data;

  if (!ad->flag || !op) {
    cd->flag = 0;
  } else {
    cd->flag = 2;
    fp_field_data_ptr p = a->field->data;
    size_t t = p->limbs;
    mp_limb_t tmp[t + 1];
    mp_limb_t qp[2];

    tmp[t] = mpn_mul_1(tmp, ad->d, t, labs(op));
    mpn_tdiv_qr(qp, cd->d, 0, tmp, t + 1, p->primelimbs, t);
    if (op < 0) {               //TODO: don't need to check c != 0 this time
      fp_neg(c, c);
    }
  }
}

static void fp_pow_mpz(element_ptr c, element_ptr a, mpz_ptr op) {
  dataptr ad = a->data;
  dataptr cd = c->data;
  if (!ad->flag) cd->flag = 0;
  else {
    mpz_t z;
    mpz_init(z);
    fp_to_mpz(z, a);
    mpz_powm(z, z, op, a->field->order);
    from_mpz(c, z);
    mpz_clear(z);
    cd->flag = 2;
  }
}

static void fp_invert(element_ptr c, element_ptr a) {
  //assumes a is invertible
  dataptr cd = c->data;
  mpz_t z;
  mpz_init(z);
  fp_to_mpz(z, a);
  mpz_invert(z, z, a->field->order);
  from_mpz(c, z);
  mpz_clear(z);
  cd->flag = 2;
}

static void fp_random(element_ptr a) {
  dataptr ad = a->data;
  mpz_t z;
  mpz_init(z);
  pbc_mpz_random(z, a->field->order);
  if (mpz_sgn(z)) {
    from_mpz(a, z);
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
  dataptr ad = a->data, bd = b->data;
  if (!ad->flag) {
    return bd->flag;
  } else {
    fp_field_data_ptr p = a->field->data;
    return mpn_cmp(ad->d, bd->d, p->limbs);
    //return memcmp(ad->d, bd->d, p->limbs);
  }
}

static int fp_sgn_odd(element_ptr a) {
  dataptr ad = a->data;
  if (!ad->flag) return 0;
  return ad->d[0] & 1 ? 1 : -1;
}

static int fp_sgn_even(element_ptr a) {
  fp_field_data_ptr p = a->field->data;
  dataptr ad = a->data;
  if (!ad->flag) return 0;
  mp_limb_t sum[p->limbs];

  int carry = mpn_add_n(sum, ad->d, ad->d, p->limbs);
  if (carry) return 1;
  return mpn_cmp(sum, p->primelimbs, p->limbs);
}


static int fp_is_sqr(element_ptr a) {
  dataptr ad = a->data;
  int res;
  mpz_t z;
  mpz_init(z);
  //0 is a square
  if (!ad->flag) return 1;
  fp_to_mpz(z, a);
  res = mpz_legendre(z, a->field->order) == 1;
  mpz_clear(z);
  return res;
}

static int fp_to_bytes(unsigned char *data, element_t a) {
  dataptr ad = a->data;
  int n = a->field->fixed_length_in_bytes;
  if (!ad->flag) {
    memset(data, 0, n);
  } else {
    mpz_t z;

    mpz_init(z);
    fp_to_mpz(z, a);
    pbc_mpz_out_raw_n(data, n, z);
    mpz_clear(z);
  }
  return n;
}

static int fp_from_bytes(element_t a, unsigned char *data) {
  dataptr ad = a->data;
  int n;
  mpz_t z;

  mpz_init(z);

  n = a->field->fixed_length_in_bytes;
  mpz_import(z, n, 1, 1, 1, 0, data);
  if (!mpz_sgn(z)) ad->flag = 0;
  else {
    ad->flag = 2;
    from_mpz(a, z);
  }
  mpz_clear(z);
  return n;
}

static void fp_out_info(FILE* str, field_ptr f) {
  element_fprintf(str, "GF(%Zd): zero flag + mpn", f->order);
}

static void fp_field_clear(field_t f) {
  fp_field_data_ptr p = f->data;
  pbc_free(p->primelimbs);
  pbc_free(p);
}

void field_init_faster_fp(field_ptr f, mpz_t prime) {
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

  f->out_info = fp_out_info;

  p = f->data = pbc_malloc(sizeof(fp_field_data_t));
  p->limbs = mpz_size(prime);
  p->bytes = p->limbs * sizeof(mp_limb_t);
  p->primelimbs = pbc_malloc(p->bytes);
  mpz_export(p->primelimbs, &p->limbs, -1, sizeof(mp_limb_t), 0, 0, prime);

  mpz_set(f->order, prime);
  f->fixed_length_in_bytes = (mpz_sizeinbase(prime, 2) + 7) / 8;
}
