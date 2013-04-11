// Quadratic extension fields.
//
// The fq_ functions are for general quadratic extensions.
// The fi_ functions are faster versions of some of these functions specialized
// for fields extended by sqrt(-1).
// TODO: Instead of lazily generating a quadratic nonresidue, in this case
// we can use sqrt(base field nqr) as the nqr of the extension.

#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdint.h> // for intptr_t
#include <stdlib.h>
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_field.h"
#include "pbc_multiz.h"
#include "pbc_fieldquadratic.h"
#include "pbc_memory.h"

// Per-element data.
typedef struct {
  // Elements have the form x + ya, where a is the square root of a quadratic
  // nonresidue in the base field.
  element_t x;
  element_t y;
} *eptr;

// Per-field data: we use ''data'' as a field_ptr to the base field.

// Return the quadratic nonresidue used to build this field.
// Should only be called from routines used exclusively by the generic quadratic
// extension code.
static inline element_ptr fq_nqr(field_ptr f) {
  return field_get_nqr((field_ptr) f->data);
}

static void fq_init(element_ptr e) {
  eptr p = e->data = pbc_malloc(sizeof(*p));
  field_ptr f = e->field->data;
  element_init(p->x, f);
  element_init(p->y, f);
}

static void fq_clear(element_ptr e) {
  eptr p = e->data;
  element_clear(p->x);
  element_clear(p->y);
  pbc_free(e->data);
}

static void fq_set_si(element_ptr e, signed long int i) {
  eptr p = e->data;
  element_set_si(p->x, i);
  element_set0(p->y);
}

static void fq_set_mpz(element_ptr e, mpz_t z) {
  eptr p = e->data;
  element_set_mpz(p->x, z);
  element_set0(p->y);
}

// Projection: attempts to convert Re(e) to mpz.
static void fq_to_mpz(mpz_t z, element_ptr e) {
  eptr p = e->data;
  element_to_mpz(z, p->x);
}

static void fq_set0(element_ptr e) {
  eptr p = e->data;
  element_set0(p->x);
  element_set0(p->y);
}

static void fq_set1(element_ptr e) {
  eptr p = e->data;
  element_set1(p->x);
  element_set0(p->y);
}

static int fq_is0(element_ptr e) {
  eptr p = e->data;
  return element_is0(p->x) && element_is0(p->y);
}

static int fq_is1(element_ptr e) {
  eptr p = e->data;
  return element_is1(p->x) && element_is0(p->y);
}

static size_t fq_out_str(FILE *stream, int base, element_ptr e) {
  size_t result = 4, status;
  eptr p = e->data;
  if (EOF == fputc('[', stream)) return 0;
  result = element_out_str(stream, base, p->x);
  if (!result) return 0;
  if (EOF == fputs(", ", stream)) return 0;
  status = element_out_str(stream, base, p->y);
  if (!status) return 0;
  if (EOF == fputc(']', stream)) return 0;
  return result + status;
}

static int fq_snprint(char *s, size_t n, element_ptr e) {
  eptr p = e->data;
  size_t result = 0, left;
  int status;

  #define clip_sub() {                   \
    result += status;                    \
    left = result >= n ? 0 : n - result; \
  }

  status = snprintf(s, n, "[");
  if (status < 0) return status;
  clip_sub();
  status = element_snprint(s + result, left, p->x);
  if (status < 0) return status;
  clip_sub();
  status = snprintf(s + result, left, ", ");
  if (status < 0) return status;
  clip_sub();
  status = element_snprint(s + result, left, p->y);
  if (status < 0) return status;
  clip_sub();
  status = snprintf(s + result, left, "]");
  if (status < 0) return status;
  return result + status;
  #undef clip_sub
}

static void fq_set_multiz(element_ptr e, multiz m) {
  eptr p = e->data;
  if (multiz_is_z(m)) {
    element_set_multiz(p->x, m);
    element_set0(p->y);
    return;
  }
  element_set_multiz(p->x, multiz_at(m, 0));
  if (2 > multiz_count(m)) element_set0(p->y);
  else element_set_multiz(p->y, multiz_at(m, 1));
}

static int fq_set_str(element_ptr e, const char *s, int base) {
  const char *cp = s;
  element_set0(e);
  while (*cp && isspace(*cp)) cp++;
  if (*cp++ != '[') return 0;
  eptr p = e->data;
  cp += element_set_str(p->x, cp, base);
  while (*cp && isspace(*cp)) cp++;
  if (*cp++ != ',') return 0;
  cp += element_set_str(p->y, cp, base);
  if (*cp++ != ']') return 0;
  return cp - s;
}

static int fq_sign(element_ptr n) {
  int res;
  eptr r = n->data;
  res = element_sign(r->x);
  if (!res) return element_sign(r->y);
  return res;
}

static void fq_add(element_ptr n, element_ptr a, element_ptr b) {
  eptr p = a->data;
  eptr q = b->data;
  eptr r = n->data;
  element_add(r->x, p->x, q->x);
  element_add(r->y, p->y, q->y);
}

static void fq_double(element_ptr n, element_ptr a) {
  eptr p = a->data;
  eptr r = n->data;
  element_double(r->x, p->x);
  element_double(r->y, p->y);
}

static void fq_sub(element_ptr n, element_ptr a, element_ptr b) {
  eptr p = a->data;
  eptr q = b->data;
  eptr r = n->data;
  element_sub(r->x, p->x, q->x);
  element_sub(r->y, p->y, q->y);
}

static void fq_set(element_ptr n, element_ptr a) {
  eptr p = a->data;
  eptr r = n->data;
  element_set(r->x, p->x);
  element_set(r->y, p->y);
}

static void fq_mul(element_ptr n, element_ptr a, element_ptr b) {
  eptr p = a->data;
  eptr q = b->data;
  eptr r = n->data;

  element_ptr nqr = fq_nqr(n->field);
  element_t e0, e1, e2;

  element_init(e0, p->x->field);
  element_init(e1, e0->field);
  element_init(e2, e0->field);
  /* naive:
  element_mul(e0, p->x, q->x);
  element_mul(e1, p->y, q->y);
  element_mul(e1, e1, nqr);
  element_add(e0, e0, e1);
  element_mul(e1, p->x, q->y);
  element_mul(e2, p->y, q->x);
  element_add(e1, e1, e2);
  element_set(r->x, e0);
  element_set(r->y, e1);
  */
  // Karatsuba:
  element_add(e0, p->x, p->y);
  element_add(e1, q->x, q->y);
  element_mul(e2, e0, e1);
  element_mul(e0, p->x, q->x);
  element_mul(e1, p->y, q->y);
  element_mul(r->x, e1, nqr);
  element_add(r->x, r->x, e0);
  element_sub(e2, e2, e0);
  element_sub(r->y, e2, e1);

  element_clear(e0);
  element_clear(e1);
  element_clear(e2);
}

static void fq_mul_mpz(element_ptr n, element_ptr a, mpz_ptr z) {
  eptr p = a->data;
  eptr r = n->data;
  element_mul_mpz(r->x, p->x, z);
  element_mul_mpz(r->y, p->y, z);
}

static void fq_mul_si(element_ptr n, element_ptr a, signed long int z) {
  eptr p = a->data;
  eptr r = n->data;
  element_mul_si(r->x, p->x, z);
  element_mul_si(r->y, p->y, z);
}

static void fq_square(element_ptr n, element_ptr a) {
  eptr p = a->data;
  eptr r = n->data;
  element_ptr nqr = fq_nqr(n->field);
  element_t e0, e1;

  element_init(e0, p->x->field);
  element_init(e1, e0->field);
  element_square(e0, p->x);
  element_square(e1, p->y);
  element_mul(e1, e1, nqr);
  element_add(e0, e0, e1);
  element_mul(e1, p->x, p->y);
  //TODO: which is faster?
  //element_add(e1, e1, e1);
  element_double(e1, e1);
  element_set(r->x, e0);
  element_set(r->y, e1);
  element_clear(e0);
  element_clear(e1);
}

static void fq_neg(element_ptr n, element_ptr a) {
  eptr p = a->data;
  eptr r = n->data;
  element_neg(r->x, p->x);
  element_neg(r->y, p->y);
}

static void fq_random(element_ptr e) {
  eptr p = e->data;
  element_random(p->x);
  element_random(p->y);
}

static int fq_cmp(element_ptr a, element_ptr b) {
  eptr p = a->data;
  eptr q = b->data;
  return element_cmp(p->x, q->x) || element_cmp(p->y, q->y);
}

static void fq_invert(element_ptr n, element_ptr a) {
  eptr p = a->data;
  eptr r = n->data;
  element_ptr nqr = fq_nqr(n->field);
  element_t e0, e1;

  element_init(e0, p->x->field);
  element_init(e1, e0->field);
  element_square(e0, p->x);
  element_square(e1, p->y);
  element_mul(e1, e1, nqr);
  element_sub(e0, e0, e1);
  element_invert(e0, e0);
  element_mul(r->x, p->x, e0);
  element_neg(e0, e0);
  element_mul(r->y, p->y, e0);

  element_clear(e0);
  element_clear(e1);
}

static void fq_from_hash(element_ptr n, void *data, int len) {
  eptr r = n->data;
  int k = len / 2;
  element_from_hash(r->x, data, k);
  element_from_hash(r->y, (char *)data + k, len - k);
}

static int fq_length_in_bytes(element_ptr e) {
  eptr p = e->data;
  return element_length_in_bytes(p->x) + element_length_in_bytes(p->y);
}

static int fq_to_bytes(unsigned char *data, element_t e) {
  eptr p = e->data;
  int len;
  len = element_to_bytes(data, p->x);
  len += element_to_bytes(data + len, p->y);
  return len;
}

static int fq_from_bytes(element_t e, unsigned char *data) {
  eptr p = e->data;
  int len;
  len = element_from_bytes(p->x, data);
  len += element_from_bytes(p->y, data + len);
  return len;
}

static int fq_is_sqr(element_ptr e) {
  //x + y sqrt(nqr) is a square iff x^2 - nqr y^2 is (in the base field)
  eptr p = e->data;
  element_t e0, e1;
  element_ptr nqr = fq_nqr(e->field);
  int result;
  element_init(e0, p->x->field);
  element_init(e1, e0->field);
  element_square(e0, p->x);
  element_square(e1, p->y);
  element_mul(e1, e1, nqr);
  element_sub(e0, e0, e1);
  result = element_is_sqr(e0);
  element_clear(e0);
  element_clear(e1);
  return result;
}

static void fq_sqrt(element_ptr n, element_ptr e) {
  eptr p = e->data;
  eptr r = n->data;
  element_ptr nqr = fq_nqr(n->field);
  element_t e0, e1, e2;

  //if (a+b sqrt(nqr))^2 = x+y sqrt(nqr) then
  //2a^2 = x +- sqrt(x^2 - nqr y^2)
  //(take the sign which allows a to exist)
  //and 2ab = y
  element_init(e0, p->x->field);
  element_init(e1, e0->field);
  element_init(e2, e0->field);
  element_square(e0, p->x);
  element_square(e1, p->y);
  element_mul(e1, e1, nqr);
  element_sub(e0, e0, e1);
  element_sqrt(e0, e0);
  //e0 = sqrt(x^2 - nqr y^2)
  element_add(e1, p->x, e0);
  element_set_si(e2, 2);
  element_invert(e2, e2);
  element_mul(e1, e1, e2);
  //e1 = (x + sqrt(x^2 - nqr y^2))/2
  if (!element_is_sqr(e1)) {
    element_sub(e1, e1, e0);
    //e1 should be a square
  }
  element_sqrt(e0, e1);
  element_add(e1, e0, e0);
  element_invert(e1, e1);
  element_mul(r->y, p->y, e1);
  element_set(r->x, e0);
  element_clear(e0);
  element_clear(e1);
  element_clear(e2);
}

static int fq_item_count(element_ptr e) {
  UNUSED_VAR(e);
  return 2;
}

static element_ptr fq_item(element_ptr e, int i) {
  eptr p = e->data;
  switch(i) {
    case 0:
      return p->x;
    case 1:
      return p->y;
    default:
      return NULL;
  }
}

static void field_clear_fq(field_ptr f) {
  UNUSED_VAR(f);
  //f->order gets cleared automatically
}

static void fq_out_info(FILE *out, field_ptr f) {
  field_ptr fbase = f->data;
  element_fprintf(out, "extension x^2 + %B, base field: ", fq_nqr(f));
  field_out_info(out, fbase);
}

// Specialized versions of some of the above for the case K[i].

static void fi_mul(element_ptr n, element_ptr a, element_ptr b) {
  eptr p = a->data;
  eptr q = b->data;
  eptr r = n->data;
  element_t e0, e1, e2;

  element_init(e0, p->x->field);
  element_init(e1, e0->field);
  element_init(e2, e0->field);
  /* Naive method:
  element_mul(e0, p->x, q->x);
  element_mul(e1, p->y, q->y);
  element_sub(e0, e0, e1);
  element_mul(e1, p->x, q->y);
  element_mul(e2, p->y, q->x);
  element_add(e1, e1, e2);
  element_set(r->x, e0);
  element_set(r->y, e1);
  */
  // Karatsuba multiplicaiton:
  element_add(e0, p->x, p->y);
  element_add(e1, q->x, q->y);
  element_mul(e2, e0, e1);
  element_mul(e0, p->x, q->x);
  element_sub(e2, e2, e0);
  element_mul(e1, p->y, q->y);
  element_sub(r->x, e0, e1);
  element_sub(r->y, e2, e1);

  element_clear(e0);
  element_clear(e1);
  element_clear(e2);
}

static void fi_square(element_ptr n, element_ptr a) {
  eptr p = a->data;
  eptr r = n->data;
  element_t e0, e1;

  element_init(e0, p->x->field);
  element_init(e1, e0->field);
  // Re(n) = x^2 - y^2 = (x+y)(x-y)
  element_add(e0, p->x, p->y);
  element_sub(e1, p->x, p->y);
  element_mul(e0, e0, e1);
  // Im(n) = 2xy
  element_mul(e1, p->x, p->y);
  element_add(e1, e1, e1);
  element_set(r->x, e0);
  element_set(r->y, e1);
  element_clear(e0);
  element_clear(e1);
}

static void fi_invert(element_ptr n, element_ptr a) {
  eptr p = a->data;
  eptr r = n->data;
  element_t e0, e1;

  element_init(e0, p->x->field);
  element_init(e1, e0->field);
  element_square(e0, p->x);
  element_square(e1, p->y);
  element_add(e0, e0, e1);
  element_invert(e0, e0);
  element_mul(r->x, p->x, e0);
  element_neg(e0, e0);
  element_mul(r->y, p->y, e0);

  element_clear(e0);
  element_clear(e1);
}

static int fi_is_sqr(element_ptr e) {
  // x + yi is a square <=> x^2 + y^2 is (in the base field).

  // Proof: (=>) if x+yi = (a+bi)^2, then a^2 - b^2 = x, 2ab = y,
  // thus (a^2 + b^2)^2 = (a^2 - b^2)^2 + (2ab)^2 =  x^2 + y^2

  // (<=) Suppose A^2 = x^2 + y^2. If there exist a, b satisfying:
  //   a^2 = (+-A + x)/2, b^2 = (+-A - x)/2
  // then (a + bi)^2 = x + yi.
  //
  // We show that exactly one of (A + x)/2, (-A + x)/2 is a quadratic residue
  // (thus a, b do exist). Suppose not. Then the product (x^2 - A^2) / 4 is
  // some quadratic residue, a contradiction since this would imply x^2 - A^2 =
  // -y^2 is also a quadratic residue, but we know -1 is not a quadratic
  // residue. QED.
  eptr p = e->data;
  element_t e0, e1;
  int result;
  element_init(e0, p->x->field);
  element_init(e1, e0->field);
  element_square(e0, p->x);
  element_square(e1, p->y);
  element_add(e0, e0, e1);
  result = element_is_sqr(e0);
  element_clear(e0);
  element_clear(e1);
  return result;
}

static void fi_sqrt(element_ptr n, element_ptr e) {
  eptr p = e->data;
  eptr r = n->data;
  element_t e0, e1, e2;

  // If (a+bi)^2 = x+yi then 2a^2 = x +- sqrt(x^2 + y^2)
  // where we choose the sign so that a exists, and 2ab = y.
  // Thus 2b^2 = - (x -+ sqrt(x^2 + y^2)).
  element_init(e0, p->x->field);
  element_init(e1, e0->field);
  element_init(e2, e0->field);
  element_square(e0, p->x);
  element_square(e1, p->y);
  element_add(e0, e0, e1);
  element_sqrt(e0, e0);
  // e0 = sqrt(x^2 + y^2)
  element_add(e1, p->x, e0);
  element_set_si(e2, 2);
  element_invert(e2, e2);
  element_mul(e1, e1, e2);
  // e1 = (x + sqrt(x^2 + y^2))/2
  if (!element_is_sqr(e1)) {
    element_sub(e1, e1, e0);
    // e1 should be a square.
  }
  element_sqrt(e0, e1);
  element_add(e1, e0, e0);
  element_invert(e1, e1);
  element_mul(r->y, p->y, e1);
  element_set(r->x, e0);
  element_clear(e0);
  element_clear(e1);
  element_clear(e2);
}

static void fi_out_info(FILE *out, field_ptr f) {
  field_ptr fbase = f->data;
  fprintf(out, "extension x^2 + 1, base field: ");
  field_out_info(out, fbase);
}

static void field_clear_fi(field_ptr f) {
  UNUSED_VAR(f);
}

// All the above should be static.

void element_field_to_quadratic(element_ptr r, element_ptr a) {
  eptr p = r->data;
  element_set(p->x, a);
  element_set0(p->y);
}

void element_field_to_fi(element_ptr a, element_ptr b) {
  element_field_to_quadratic(a, b);
}

static element_ptr fq_get_x(element_ptr a) {
  return ((eptr) a->data)->x;
}

static element_ptr fq_get_y(element_ptr a) {
  return ((eptr) a->data)->y;
}

void field_init_quadratic(field_ptr f, field_ptr fbase) {
  field_init(f);

  f->field_clear = field_clear_fq;
  f->data = fbase;

  f->init = fq_init;
  f->clear = fq_clear;
  f->set_si = fq_set_si;
  f->set_mpz = fq_set_mpz;
  f->to_mpz = fq_to_mpz;
  f->out_str = fq_out_str;
  f->snprint = fq_snprint;
  f->set_multiz = fq_set_multiz;
  f->set_str = fq_set_str;
  f->sign = fq_sign;
  f->add = fq_add;
  f->sub = fq_sub;
  f->set = fq_set;
  f->mul = fq_mul;
  f->mul_mpz = fq_mul_mpz;
  f->mul_si = fq_mul_si;
  f->square = fq_square;
  f->doub = fq_double;
  f->neg = fq_neg;
  f->cmp = fq_cmp;
  f->invert = fq_invert;
  f->random = fq_random;
  f->from_hash = fq_from_hash;
  f->is1 = fq_is1;
  f->is0 = fq_is0;
  f->set0 = fq_set0;
  f->set1 = fq_set1;
  f->is_sqr = fq_is_sqr;
  f->sqrt = fq_sqrt;
  f->to_bytes = fq_to_bytes;
  f->from_bytes = fq_from_bytes;
  f->out_info = fq_out_info;
  f->item_count = fq_item_count;
  f->item = fq_item;
  f->get_x = fq_get_x;
  f->get_y = fq_get_y;

  mpz_mul(f->order, fbase->order, fbase->order);
  if (fbase->fixed_length_in_bytes < 0) {
    f->length_in_bytes = fq_length_in_bytes;
    f->fixed_length_in_bytes = -1;
  } else {
    f->fixed_length_in_bytes = 2 * fbase->fixed_length_in_bytes;
  }
}

void field_init_fi(field_ptr f, field_ptr fbase) {
  field_init(f);
  f->field_clear = field_clear_fi;
  f->data = fbase;
  f->init = fq_init;
  f->clear = fq_clear;
  f->set_si = fq_set_si;
  f->set_mpz = fq_set_mpz;
  f->to_mpz = fq_to_mpz;
  f->out_str = fq_out_str;
  f->snprint = fq_snprint;
  f->set_multiz = fq_set_multiz;
  f->set_str = fq_set_str;
  f->sign = fq_sign;
  f->add = fq_add;
  f->sub = fq_sub;
  f->set = fq_set;
  f->mul = fi_mul;
  f->mul_mpz = fq_mul_mpz;
  f->mul_si = fq_mul_si;
  f->square = fi_square;
  f->doub = fq_double;
  f->neg = fq_neg;
  f->cmp = fq_cmp;
  f->invert = fi_invert;
  f->random = fq_random;
  f->from_hash = fq_from_hash;
  f->is1 = fq_is1;
  f->is0 = fq_is0;
  f->set0 = fq_set0;
  f->set1 = fq_set1;
  f->is_sqr = fi_is_sqr;
  f->sqrt = fi_sqrt;
  f->to_bytes = fq_to_bytes;
  f->from_bytes = fq_from_bytes;
  f->out_info = fi_out_info;
  f->item_count = fq_item_count;
  f->item = fq_item;
  f->get_x = fq_get_x;
  f->get_y = fq_get_y;

  mpz_mul(f->order, fbase->order, fbase->order);
  if (fbase->fixed_length_in_bytes < 0) {
    f->length_in_bytes = fq_length_in_bytes;
    f->fixed_length_in_bytes = -1;
  } else {
    f->fixed_length_in_bytes = 2 * fbase->fixed_length_in_bytes;
  }
}
