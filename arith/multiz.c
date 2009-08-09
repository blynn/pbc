// Multinomials over Z.
// e.g. [[1, 2], 3, [4, [5, 6]]] means
// (1 + 2y) + 3 x + (4 + (5 + 6z)y)x^2
// Convenient interchange format for different groups, rings, and fields.

// TODO: Canonicalize, e.g. [[1]], 0, 0] --> 1.

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_field.h"
#include "pbc_multiz.h"
#include "pbc_random.h"
#include "pbc_fp.h"
#include "pbc_memory.h"
#include "misc/darray.h"

// Per-element data.
typedef struct {
  // Either it's an mpz, or a list of mpzs.
  char type;
  union {
    mpz_t z;
    darray_t a;
  };
} *multiz;

enum {
  T_MPZ,
  T_ARR,
};

static multiz multiz_new_empty_list(void) {
  multiz ep = pbc_malloc(sizeof(*ep));
  ep->type = T_ARR;
  darray_init(ep->a);
  return ep;
}

void multiz_append(element_ptr x, element_ptr e) {
  multiz l = x->data;
  darray_append(l->a, e->data);
}

static multiz multiz_new(void) {
  multiz ep = pbc_malloc(sizeof(*ep));
  ep->type = T_MPZ;
  mpz_init(ep->z);
  return ep;
}

static void f_init(element_ptr e) {
  e->data = multiz_new();
}

static void multiz_free(multiz ep) {
  switch(ep->type) {
    case T_MPZ:
      mpz_clear(ep->z);
      break;
    default:
      PBC_ASSERT(T_ARR == ep->type, "no such type");
      void clearit(void *data) { multiz_free(data); }
      darray_forall(ep->a, clearit);
      darray_clear(ep->a);
      break;
  }
  pbc_free(ep);
}

static void f_clear(element_ptr e) {
  multiz_free(e->data);
}

element_ptr multiz_new_list(element_ptr e) {
  element_ptr x = pbc_malloc(sizeof(*x));
  element_init_same_as(x, e);
  multiz_free(x->data);
  x->data = multiz_new_empty_list();
  multiz_append(x, e);
  return x;
}

static void f_set_si(element_ptr e, signed long int op) {
  multiz_free(e->data);
  f_init(e);
  multiz ep = e->data;
  mpz_set_si(ep->z, op);
}

static void f_set_mpz(element_ptr e, mpz_ptr z) {
  multiz_free(e->data);
  f_init(e);
  multiz ep = e->data;
  mpz_set(ep->z, z);
}

static void f_set0(element_ptr e) {
  multiz_free(e->data);
  f_init(e);
}

static void f_set1(element_ptr e) {
  multiz_free(e->data);
  f_init(e);
  multiz ep = e->data;
  mpz_set_ui(ep->z, 1);
}

static size_t multiz_out_str(FILE *stream, int base, multiz ep) {
  switch(ep->type) {
    case T_MPZ:
      return mpz_out_str(stream, base, ep->z);
    default:
      PBC_ASSERT(T_ARR == ep->type, "no such type");
      fputc('[', stream);
      size_t res = 1;
      int n = darray_count(ep->a);
      int i;
      for(i = 0; i < n; i++) {
	if (i) res += 2, fputs(", ", stream);
	res += multiz_out_str(stream, base, darray_at(ep->a, i));
      }
      fputc(']', stream);
      res++;
      return res;
  }
}

static size_t f_out_str(FILE *stream, int base, element_ptr e) {
  return multiz_out_str(stream, base, e->data);
}

void multiz_to_mpz(mpz_ptr z, multiz ep) {
  while(ep->type == T_ARR) ep = darray_at(ep->a, 0);
  PBC_ASSERT(T_MPZ == ep->type, "no such type");
  mpz_set(z, ep->z);
}

static void f_to_mpz(mpz_ptr z, element_ptr a) {
  multiz_to_mpz(z, a->data);
}

static int multiz_sgn(multiz ep) {
  while(ep->type == T_ARR) ep = darray_at(ep->a, 0);
  PBC_ASSERT(T_MPZ == ep->type, "no such type");
  return mpz_sgn(ep->z);
}

static int f_sgn(element_ptr a) {
  return multiz_sgn(a->data);
}

static multiz multiz_unary(const multiz y,
    void (*fun)(mpz_ptr, const mpz_ptr)) {
  multiz x = pbc_malloc(sizeof(*x));
  switch(y->type) {
    case T_MPZ:
      x->type = T_MPZ;
      mpz_init(x->z);
      fun(x->z, y->z);
      break;
    default:
      PBC_ASSERT(T_ARR == ep->type, "no such type");
      x->type = T_ARR;
      darray_init(x->a);
      void add_to_x(void *data) {
	darray_append(x->a, multiz_unary(data, fun));
      }
      darray_forall(y->a, add_to_x);
      break;
  }
  return x;
}

// Need this wrapper to suppress a warning; I can't figure out what
// mpz_set()'s type is. Same for other our_* functions.
static void our_mpz_set(mpz_ptr x, const mpz_ptr a) { mpz_set(x, a); }
static multiz multiz_clone(multiz y) {
  return multiz_unary(y, our_mpz_set);
}

static multiz multiz_new_add(multiz a, multiz b) {
  if (T_MPZ == a->type) {
    if (T_MPZ == b->type) {
      multiz x = multiz_new();
      mpz_add(x->z, a->z, b->z);
      return x;
    } else {
      multiz x = multiz_clone(b);
      multiz z = x;
      PBC_ASSERT(T_ARR == z->type, "no such type");
      while(z->type == T_ARR) z = darray_at(z->a, 0);
      mpz_add(z->z, a->z, z->z);
      return x;
    }
  } else {
    PBC_ASSERT(T_ARR == a->type, "no such type");
    if (T_MPZ == b->type) {
      multiz x = multiz_clone(a);
      multiz z = x;
      PBC_ASSERT(T_ARR == z->type, "no such type");
      while(z->type == T_ARR) z = darray_at(z->a, 0);
      mpz_add(z->z, b->z, z->z);
      return x;
    } else {
      PBC_ASSERT(T_ARR == b->type, "no such type");
      int m = darray_count(a->a);
      int n = darray_count(b->a);
      int min = m < n ? m : n;
      int max = m > n ? m : n;
      multiz x = multiz_new_empty_list();
      int i;
      for(i = 0; i < min; i++) {
	multiz z = multiz_new_add(darray_at(a->a, i), darray_at(b->a, i));
	darray_append(x->a, z);
      }
      for(; i < max; i++) {
	darray_append(x->a, multiz_clone(darray_at(m > n ? a->a : b->a, i)));
      }
      return x;
    }
  }
}

static void f_add(element_ptr n, element_ptr a, element_ptr b) {
  void *delme = n->data;
  n->data = multiz_new_add(a->data, b->data);
  multiz_free(delme);
}

static void z_sub(element_ptr n, element_ptr a, element_ptr b) {
  mpz_sub(n->data, a->data, b->data);
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

static void our_mpz_neg(mpz_ptr x, const mpz_ptr a) { mpz_neg(x, a); }
static multiz multiz_new_neg(multiz z) {
  return multiz_unary(z, our_mpz_neg);
}

static void f_set(element_ptr n, element_ptr a) {
  void *delme = n->data;
  n->data = multiz_clone(a->data);
  multiz_free(delme);
}

static void f_neg(element_ptr n, element_ptr a) {
  void *delme = n->data;
  n->data = multiz_new_neg(a->data);
  multiz_free(delme);
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

static void f_field_clear(field_t f) {
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

static int z_length_in_bytes(element_ptr a) {
  return (mpz_sizeinbase(a->data, 2) + 7) / 8 + 4;
}

static void z_out_info(FILE *out, field_ptr f) {
  UNUSED_VAR(f);
  fprintf(out, "Z: wrapped GMP");
}

static int f_set_str(element_ptr e, const char *s, int base) {
  // TODO: Square brackets.
  mpz_t z;
  mpz_init(z);
  int result = pbc_mpz_set_str(z, s, base);
  f_set_mpz(e, z);
  mpz_clear(z);
  return result;
}

static void f_set_multiz(element_ptr e, element_ptr m) {
  multiz_free(e->data);
  e->data = multiz_clone(m->data);;
}

void field_init_multiz(field_ptr f) {
  field_init(f);
  f->init = f_init;
  f->clear = f_clear;
  f->set_si = f_set_si;
  f->set_mpz = f_set_mpz;
  f->set_multiz = f_set_multiz;
  f->set_str = f_set_str;
  f->out_str = f_out_str;
  f->sign = f_sgn;
  f->add = f_add;
  f->sub = z_sub;
  f->set = f_set;
  f->mul = z_mul;
  f->mul_mpz = z_mul_mpz;
  f->mul_si = z_mul_si;
  f->pow_mpz = z_pow_mpz;
  f->neg = f_neg;
  f->cmp = z_cmp;
  f->invert = z_invert;
  f->div = z_div;
  f->random = z_random;
  f->from_hash = z_from_hash;
  f->is1 = z_is1;
  f->is0 = z_is0;
  f->set0 = f_set0;
  f->set1 = f_set1;
  f->field_clear = f_field_clear;
  f->to_bytes = z_to_bytes;
  f->from_bytes = z_from_bytes;
  f->to_mpz = f_to_mpz;
  f->length_in_bytes = z_length_in_bytes;

  f->out_info = z_out_info;

  mpz_set_ui(f->order, 0);
  f->data = NULL;
  f->fixed_length_in_bytes = -1;
}
