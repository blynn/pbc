// Multinomials over Z.
// e.g. [[1, 2], 3, [4, [5, 6]]] means
// (1 + 2y) + 3 x + (4 + (5 + 6z)y)x^2
// Convenient interchange format for different groups, rings, and fields.

// TODO: Canonicalize, e.g. [[1]], 0, 0] --> 1.

#include <stdarg.h>
#include <stdio.h>
#include <stdint.h> // for intptr_t
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
struct multiz_s {
  // Either it's an mpz, or a list of mpzs.
  char type;
  union {
    mpz_t z;
    darray_t a;
  };
};

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
      darray_forall(ep->a, (void(*)(void*))multiz_free);
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

static void add_to_x(void *data,
                     multiz x,
                     void (*fun)(mpz_t, const mpz_t, void *scope_ptr),
                     void *scope_ptr);

static multiz multiz_new_unary(const multiz y,
    void (*fun)(mpz_t, const mpz_t, void *scope_ptr), void *scope_ptr) {
  multiz x = pbc_malloc(sizeof(*x));
  switch(y->type) {
    case T_MPZ:
      x->type = T_MPZ;
      mpz_init(x->z);
      fun(x->z, y->z, scope_ptr);
      break;
    default:
      PBC_ASSERT(T_ARR == ep->type, "no such type");
      x->type = T_ARR;
      darray_init(x->a);
      darray_forall4(y->a,
                     (void(*)(void*,void*,void*,void*))add_to_x,
                     x,
                     fun,
                     scope_ptr);
      break;
  }
  return x;
}

static void add_to_x(void *data,
                     multiz x,
                     void (*fun)(mpz_t, const mpz_t, void *scope_ptr),
                     void *scope_ptr) {
  darray_append(x->a, multiz_new_unary(data, fun, scope_ptr));
}

static void mpzset(mpz_t dst, const mpz_t src, void *scope_ptr) {
  UNUSED_VAR(scope_ptr);
  mpz_set(dst, src);
}

static multiz multiz_clone(multiz y) {
  return multiz_new_unary(y, (void(*)(mpz_t, const mpz_t, void *))mpzset, NULL);
}

static multiz multiz_new_bin(const multiz a, const multiz b,
    void (*fun)(mpz_t, const mpz_t, const mpz_t)) {
  if (T_MPZ == a->type) {
    if (T_MPZ == b->type) {
      multiz x = multiz_new();
      fun(x->z, a->z, b->z);
      return x;
    } else {
      multiz x = multiz_clone(b);
      multiz z = x;
      PBC_ASSERT(T_ARR == z->type, "no such type");
      while(z->type == T_ARR) z = darray_at(z->a, 0);
      fun(z->z, a->z, z->z);
      return x;
    }
  } else {
    PBC_ASSERT(T_ARR == a->type, "no such type");
    if (T_MPZ == b->type) {
      multiz x = multiz_clone(a);
      multiz z = x;
      PBC_ASSERT(T_ARR == z->type, "no such type");
      while(z->type == T_ARR) z = darray_at(z->a, 0);
      fun(z->z, b->z, z->z);
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
        multiz z = multiz_new_bin(darray_at(a->a, i), darray_at(b->a, i), fun);
        darray_append(x->a, z);
      }
      multiz zero = multiz_new();
      for(; i < max; i++) {
        multiz z = multiz_new_bin(m > n ? darray_at(a->a, i) : zero,
                                  n > m ? darray_at(b->a, i) : zero,
                                  fun);
        darray_append(x->a, z);
      }
      multiz_free(zero);
      return x;
    }
  }
}
static multiz multiz_new_add(const multiz a, const multiz b) {
  return multiz_new_bin(a, b, mpz_add);
}

static void f_add(element_ptr n, element_ptr a, element_ptr b) {
  multiz delme = n->data;
  n->data = multiz_new_add(a->data, b->data);
  multiz_free(delme);
}

static multiz multiz_new_sub(const multiz a, const multiz b) {
  return multiz_new_bin(a, b, mpz_sub);
}
static void f_sub(element_ptr n, element_ptr a, element_ptr b) {
  multiz delme = n->data;
  n->data = multiz_new_sub(a->data, b->data);
  multiz_free(delme);
}

static void mpzmul(mpz_t x, const mpz_t y, const mpz_t z) {
  mpz_mul(x, y, z);
}

static multiz multiz_new_mul(const multiz a, const multiz b) {
  if (T_MPZ == a->type) {
    // Multiply each coefficient of b by a->z.
    return multiz_new_unary(b, (void(*)(mpz_t, const mpz_t, void *))mpzmul, a->z);
  } else {
    PBC_ASSERT(T_ARR == a->type, "no such type");
    if (T_MPZ == b->type) {
      // Multiply each coefficient of a by b->z.
      return multiz_new_unary(a, (void(*)(mpz_t, const mpz_t, void *))mpzmul, b->z);
    } else {
      PBC_ASSERT(T_ARR == b->type, "no such type");
      int m = darray_count(a->a);
      int n = darray_count(b->a);
      int max = m + n - 1;
      multiz x = multiz_new_empty_list();
      int i;
      multiz zero = multiz_new();
      for(i = 0; i < max; i++) {
        multiz z = multiz_new();
        int j;
        for (j = 0; j <= i; j++) {
          multiz y = multiz_new_mul(j < m ? darray_at(a->a, j) : zero,
                                    i - j < n ? darray_at(b->a, i - j) : zero);
          multiz t = multiz_new_add(z, y);
          multiz_free(y);
          multiz_free(z);
          z = t;
        }
        darray_append(x->a, z);
      }
      multiz_free(zero);
      return x;
    }
  }
}
static void f_mul(element_ptr n, element_ptr a, element_ptr b) {
  multiz delme = n->data;
  n->data = multiz_new_mul(a->data, b->data);
  multiz_free(delme);
}

static void f_mul_mpz(element_ptr n, element_ptr a, mpz_ptr z) {
  multiz delme = n->data;
  n->data = multiz_new_unary(a->data, (void(*)(mpz_t, const mpz_t, void *))mpzmul, z);
  multiz_free(delme);
}

static void mulsi(mpz_t x, const mpz_t y, signed long *i) {
  mpz_mul_si(x, y, *i);
}

static void f_mul_si(element_ptr n, element_ptr a, signed long int z) {
  multiz delme = n->data;
  n->data = multiz_new_unary(a->data, (void(*)(mpz_t, const mpz_t, void *))mulsi, &z);
  multiz_free(delme);
}

static void mpzneg(mpz_t dst, const mpz_t src, void *scope_ptr) {
  UNUSED_VAR(scope_ptr);
  mpz_neg(dst, src);
}

static multiz multiz_new_neg(multiz z) {
  return multiz_new_unary(z, (void(*)(mpz_t, const mpz_t, void *))mpzneg, NULL);
}

static void f_set(element_ptr n, element_ptr a) {
  multiz delme = n->data;
  n->data = multiz_clone(a->data);
  multiz_free(delme);
}

static void f_neg(element_ptr n, element_ptr a) {
  multiz delme = n->data;
  n->data = multiz_new_neg(a->data);
  multiz_free(delme);
}

static void f_div(element_ptr c, element_ptr a, element_ptr b) {
  mpz_t d;
  mpz_init(d);
  element_to_mpz(d, b);
  multiz delme = c->data;
  c->data = multiz_new_unary(a->data, (void(*)(mpz_t, const mpz_t, void *))mpz_tdiv_q, d);
  mpz_clear(d);
  multiz_free(delme);
}

// Doesn't make sense if order is infinite.
static void f_random(element_ptr n) {
  multiz delme = n->data;
  f_init(n);
  multiz_free(delme);
}

static void f_from_hash(element_ptr n, void *data, int len) {
  mpz_t z;
  mpz_init(z);
  mpz_import(z, len, -1, 1, -1, 0, data);
  f_set_mpz(n, z);
  mpz_clear(z);
}

static int f_is1(element_ptr n) {
  multiz ep = n->data;
  return ep->type == T_MPZ && !mpz_cmp_ui(ep->z, 1);
}

int multiz_is0(multiz m) {
  return m->type == T_MPZ && mpz_is0(m->z);
}

static int f_is0(element_ptr n) {
  return multiz_is0(n->data);
}

static int f_item_count(element_ptr e) {
  multiz z = e->data;
  if (T_MPZ == z->type) return 0;
  return darray_count(z->a);
}

// TODO: Redesign multiz so this doesn't leak.
static element_ptr f_item(element_ptr e, int i) {
  multiz z = e->data;
  if (T_MPZ == z->type) return NULL;
  element_ptr r = malloc(sizeof(*r));
  r->field = e->field;
  r->data = darray_at(z->a, i);
  return r;
}

// Usual meaning when both are integers.
// Otherwise, compare coefficients.
static int multiz_cmp(multiz a, multiz b) {
  if (T_MPZ == a->type) {
    if (T_MPZ == b->type) {
      // Simplest case: both are integers.
      return mpz_cmp(a->z, b->z);
    }
    // Leading coefficient of b.
    while(T_ARR == b->type) b = darray_last(b->a);
    PBC_ASSERT(T_MPZ == b->type, "no such type");
    return -mpz_sgn(b->z);
  }
  PBC_ASSERT(T_ARR == a->type, "no such type");
  if (T_MPZ == b->type) {
    // Leading coefficient of a.
    while(T_ARR == a->type) a = darray_last(a->a);
    PBC_ASSERT(T_MPZ == a->type, "no such type");
    return mpz_sgn(a->z);
  }
  PBC_ASSERT(T_ARR == b->type, "no such type");
  int m = darray_count(a->a);
  int n = darray_count(b->a);
  if (m > n) {
    // Leading coefficient of a.
    while(T_ARR == a->type) a = darray_last(a->a);
    PBC_ASSERT(T_MPZ == a->type, "no such type");
    return mpz_sgn(a->z);
  }
  if (n > m) {
    // Leading coefficient of b.
    while(T_ARR == b->type) b = darray_last(b->a);
    PBC_ASSERT(T_MPZ == b->type, "no such type");
    return -mpz_sgn(b->z);
  }
  for(n--; n >= 0; n--) {
    int i = multiz_cmp(darray_at(a->a, n), darray_at(b->a, n));
    if (i) return i;
  }
  return 0;
}
static int f_cmp(element_ptr x, element_ptr y) {
  return multiz_cmp(x->data, y->data);
}

static void f_field_clear(field_t f) { UNUSED_VAR (f); }

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

static void f_out_info(FILE *out, field_ptr f) {
  UNUSED_VAR(f);
  fprintf(out, "Z multinomials");
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

static void f_set_multiz(element_ptr e, multiz m) {
  multiz delme = e->data;
  e->data = multiz_clone(m);
  multiz_free(delme);
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
  f->sub = f_sub;
  f->set = f_set;
  f->mul = f_mul;
  f->mul_mpz = f_mul_mpz;
  f->mul_si = f_mul_si;
  f->neg = f_neg;
  f->cmp = f_cmp;
  f->div = f_div;
  f->random = f_random;
  f->from_hash = f_from_hash;
  f->is1 = f_is1;
  f->is0 = f_is0;
  f->set0 = f_set0;
  f->set1 = f_set1;
  f->field_clear = f_field_clear;
  f->to_bytes = z_to_bytes;
  f->from_bytes = z_from_bytes;
  f->to_mpz = f_to_mpz;
  f->length_in_bytes = z_length_in_bytes;
  f->item = f_item;
  f->item_count = f_item_count;

  f->out_info = f_out_info;

  mpz_set_ui(f->order, 0);
  f->data = NULL;
  f->fixed_length_in_bytes = -1;
}

int multiz_is_z(multiz m) {
  return T_MPZ == m->type;
}

int multiz_count(multiz m) {
  if (T_ARR != m->type) return -1;
  return darray_count(m->a);
}

multiz multiz_at(multiz m, int i) {
  PBC_ASSERT(T_ARR == m->type, "wrong type");
  PBC_ASSERT(darray_count(m->a) > i, "out of bounds");
  return darray_at(m->a, i);
}
