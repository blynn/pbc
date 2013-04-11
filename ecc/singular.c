#include <stdarg.h>
#include <stdio.h>
#include <stdint.h> // for intptr_t
#include <stdlib.h>
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_field.h"
#include "pbc_curve.h"
#include "pbc_param.h"
#include "pbc_pairing.h"
#include "pbc_fp.h"
#include "pbc_memory.h"

//TODO: Store as integer mod ring instead and convert at last minute?
struct point_s {
  int inf_flag;
  element_t x;
  element_t y;
};
typedef struct point_s *point_ptr;
typedef struct point_s point_t[1];

static void sn_init(element_ptr e) {
  field_ptr f = e->field->data;
  e->data = pbc_malloc(sizeof(point_t));
  point_ptr p = e->data;
  element_init(p->x, f);
  element_init(p->y, f);
  p->inf_flag = 1;
}

static void sn_clear(element_ptr e) {
  point_ptr p = e->data;
  element_clear(p->x);
  element_clear(p->y);
  pbc_free(e->data);
}

static void sn_set0(element_ptr x) {
  point_ptr p = x->data;
  p->inf_flag = 1;
}

static int sn_is0(element_ptr x) {
  point_ptr p = x->data;
  return p->inf_flag;
}

//singular with node: y^2 = x^3 + x^2
static void sn_random(element_t a) {
  point_ptr p = a->data;
  element_t t;

  element_init(t, p->x->field);
  p->inf_flag = 0;
  do {
    element_random(p->x);
    if (element_is0(p->x)) continue;
    element_square(t, p->x);
    element_add(t, t, p->x);
    element_mul(t, t, p->x);
  } while (!element_is_sqr(t));
  element_sqrt(p->y, t);

  element_clear(t);
}

static inline void sn_double_no_check(point_ptr r, point_ptr p) {
  element_t lambda, e0, e1;

  element_init(lambda, p->x->field);
  element_init(e0, p->x->field);
  element_init(e1, p->x->field);
  //same point: double them

  //lambda = (3x^2 + 2x) / 2y
  element_mul_si(lambda, p->x, 3);
  element_set_si(e0, 2);
  element_add(lambda, lambda, e0);
  element_mul(lambda, lambda, p->x);
  element_add(e0, p->y, p->y);
  element_invert(e0, e0);
  element_mul(lambda, lambda, e0);
  //x1 = lambda^2 - 2x - 1
  element_add(e1, p->x, p->x);
  element_square(e0, lambda);
  element_sub(e0, e0, e1);
  element_set_si(e1, 1);
  element_sub(e0, e0, e1);
  //y1 = (x - x1)lambda - y
  element_sub(e1, p->x, e0);
  element_mul(e1, e1, lambda);
  element_sub(e1, e1, p->y);

  element_set(r->x, e0);
  element_set(r->y, e1);
  r->inf_flag = 0;

  element_clear(lambda);
  element_clear(e0);
  element_clear(e1);
  return;
}

static void sn_double(element_t c, element_t a) {
  point_ptr r = c->data;
  point_ptr p = a->data;
  if (p->inf_flag) {
    r->inf_flag = 1;
    return;
  }
  if (element_is0(p->y)) {
    r->inf_flag = 1;
    return;
  }
  sn_double_no_check(r, p);
}

static void sn_set(element_ptr c, element_ptr a) {
  point_ptr r = c->data, p = a->data;
  if (p->inf_flag) {
    r->inf_flag = 1;
    return;
  }
  r->inf_flag = 0;
  element_set(r->x, p->x);
  element_set(r->y, p->y);
}

static void sn_add(element_t c, element_t a, element_t b) {
  point_ptr r = c->data;
  point_ptr p = a->data;
  point_ptr q = b->data;
  if (p->inf_flag) {
    sn_set(c, b);
    return;
  }
  if (q->inf_flag) {
    sn_set(c, a);
    return;
  }
  if (!element_cmp(p->x, q->x)) {
    if (!element_cmp(p->y, q->y)) {
      if (element_is0(p->y)) {
        r->inf_flag = 1;
        return;
      } else {
        sn_double_no_check(r, p);
        return;
      }
    }
    //points are inverses of each other
    r->inf_flag = 1;
    return;
  } else {
    element_t lambda, e0, e1;

    element_init(lambda, p->x->field);
    element_init(e0, p->x->field);
    element_init(e1, p->x->field);

    //lambda = (y2-y1)/(x2-x1)
    element_sub(e0, q->x, p->x);
    element_invert(e0, e0);
    element_sub(lambda, q->y, p->y);
    element_mul(lambda, lambda, e0);
    //x3 = lambda^2 - x1 - x2 - 1
    element_square(e0, lambda);
    element_sub(e0, e0, p->x);
    element_sub(e0, e0, q->x);
    element_set1(e1);
    element_sub(e0, e0, e1);
    //y3 = (x1-x3)lambda - y1
    element_sub(e1, p->x, e0);
    element_mul(e1, e1, lambda);
    element_sub(e1, e1, p->y);

    element_set(r->x, e0);
    element_set(r->y, e1);
    r->inf_flag = 0;

    element_clear(lambda);
    element_clear(e0);
    element_clear(e1);
  }
}

static void sn_invert(element_ptr c, element_ptr a) {
  point_ptr r = c->data, p = a->data;

  if (p->inf_flag) {
    r->inf_flag = 1;
    return;
  }
  r->inf_flag = 0;
  element_set(r->x, p->x);
  element_neg(r->y, p->y);
}

static void sn_field_clear(field_ptr c) {
  UNUSED_VAR(c);
}

/* TODO: Write a test program that uses these functions.

// Nonsingular points on sn curves map to finite field elements via
//   (x, y) --> (y + x)/(y - x)
// The reverse map is
//   a --> (4a/(a-1)^2, 4a(a+1)/(a-1)^3)

void sn_point_to_field(element_t out, point_ptr P) {
  element_t e0, e1;
  if (P->inf_flag) {
    element_set1(out);
    return;
  }
  element_init(e0, out->field);
  element_init(e1, out->field);
  element_add(e0, P->y, P->x);
  element_sub(e1, P->y, P->x);
  element_invert(e1, e1);
  element_mul(out, e0, e1);
  element_clear(e0);
  element_clear(e1);
}

static void sn_field_to_point(point_ptr P, element_t in) {
  element_t e0, e1, e2;

  if (element_is1(in)) {
    P->inf_flag = 1;
    return;
  }
  element_init(e0, in->field);
  element_init(e1, in->field);
  element_init(e2, in->field);

  element_set1(e1);
  element_sub(e0, in, e1);
  element_invert(e0, e0);

  element_mul_si(e2, in, 4);

  element_add(P->y, in, e1);

  element_mul(e1, e0, e0);
  element_mul(P->x, e1, e2);
  element_mul(P->y, P->y, e2);
  element_mul(P->y, P->y, e0);
  element_mul(P->y, P->y, e1);
  P->inf_flag = 0;

  element_clear(e0);
  element_clear(e1);
  element_clear(e2);
}
*/

static size_t sn_out_str(FILE *stream, int base, element_ptr a) {
  point_ptr p = a->data;
  size_t result, status;
  if (p->inf_flag) {
    if (EOF == fputc('O', stream)) return 0;
    return 1;
  }
  result = element_out_str(stream, base, p->x);
  if (!result) return 0;
  if (EOF == fputc(' ', stream)) return 0;
  status = element_out_str(stream, base, p->y);
  if (!status) return 0;
  return result + status + 1;
}

void naive_generic_pow_mpz(element_ptr x, element_ptr a, mpz_ptr n);
void field_init_curve_singular_with_node(field_t c, field_t field) {
  mpz_set(c->order, field->order);
  c->data = (void *) field;
  c->init = sn_init;
  c->clear = sn_clear;
  c->random = sn_random;
  //c->from_x = cc_from_x;
  //c->from_hash = cc_from_hash;
  c->set = sn_set;
  c->invert = c->neg = sn_invert;
  c->square = c->doub = sn_double;
  c->mul = c->add = sn_add;
  c->set1 = c->set0 = sn_set0;
  c->is1 = c->is0 = sn_is0;
  c->mul_mpz = element_pow_mpz;
  c->out_str = sn_out_str;
  c->field_clear = sn_field_clear;
}

//TODO: the following code is useless as the Tate pairing is degenerate on singular curves
static void sn_miller(element_t res, mpz_t q, element_t P,
    element_ptr Qx, element_ptr Qy) {
  //collate divisions
  int m;
  element_t v, vd;
  element_t Z;
  element_t a, b, c;
  element_t e0, e1;
  element_ptr Zx;
  element_ptr Zy;
  const element_ptr Px = curve_x_coord(P);
  const element_ptr Py = curve_y_coord(P);

  #define do_vertical(e)     \
    element_sub(e0, Qx, Zx); \
    element_mul(e, e, e0);

  //a = -slope_tangent(Z.x, Z.y);
  //b = 1;
  //c = -(Z.y + a * Z.x);
  //but we multiply by 2*Z.y to avoid division
  //a = -Zx * (Zx + Zx + Zx + 2)
  //b = 2 * Zy
  //c = -(2 Zy^2 + a Zx);
  #define do_tangent(e)      \
    element_double(e0, Zx);  \
    element_add(a, Zx, e0);  \
    element_set_si(e0, 2);   \
    element_add(a, a, e0);   \
    element_mul(a, a, Zx);   \
    element_neg(a, a);       \
    element_add(b, Zy, Zy);  \
    element_mul(e0, b, Zy);  \
    element_mul(c, a, Zx);   \
    element_add(c, c, e0);   \
    element_neg(c, c);       \
    element_mul(e0, a, Qx);  \
    element_mul(e1, b, Qy);  \
    element_add(e0, e0, e1); \
    element_add(e0, e0, c);  \
    element_mul(e, e, e0);

  //a = -(B.y - A.y) / (B.x - A.x);
  //b = 1;
  //c = -(A.y + a * A.x);
  //but we'll multiply by B.x - A.x to avoid division
  #define do_line(e)         \
    element_sub(b, Px, Zx);  \
    element_sub(a, Zy, Py);  \
    element_mul(e0, b, Zy);  \
    element_mul(c, a, Zx);   \
    element_add(c, c, e0);   \
    element_neg(c, c);       \
    element_mul(e0, a, Qx);  \
    element_mul(e1, b, Qy);  \
    element_add(e0, e0, e1); \
    element_add(e0, e0, c);  \
    element_mul(e, e, e0);

  element_init(a, Px->field);
  element_init(b, Px->field);
  element_init(c, Px->field);
  element_init(e0, res->field);
  element_init(e1, res->field);

  element_init(v, res->field);
  element_init(vd, res->field);
  element_init(Z, P->field);

  element_set(Z, P);
  Zx = curve_x_coord(Z);
  Zy = curve_y_coord(Z);

  element_set1(v);
  element_set1(vd);
  m = mpz_sizeinbase(q, 2) - 2;

  while(m >= 0) {
    element_mul(v, v, v);
    element_mul(vd, vd, vd);
    do_tangent(v);
    element_double(Z, Z);
    do_vertical(vd);
    if (mpz_tstbit(q, m)) {
      do_line(v);
      element_add(Z, Z, P);
      do_vertical(vd);
    }
    m--;
  }
  #undef do_tangent
  #undef do_vertical
  #undef do_line

  element_invert(vd, vd);
  element_mul(res, v, vd);

  element_clear(v);
  element_clear(vd);
  element_clear(Z);
  element_clear(a);
  element_clear(b);
  element_clear(c);
  element_clear(e0);
  element_clear(e1);
}

struct sn_pairing_data_s {
  field_t Fq, Eq;
};
typedef struct sn_pairing_data_s sn_pairing_data_t[1];
typedef struct sn_pairing_data_s *sn_pairing_data_ptr;

static void sn_pairing(element_ptr out, element_ptr in1, element_ptr in2,
    pairing_t pairing) {
  sn_pairing_data_ptr p = pairing->data;
  element_ptr Q = in2;
  element_t e0;
  element_t R, QR;
  element_init(R, p->Eq);
  element_init(QR, p->Eq);
  element_random(R);
  element_init(e0, out->field);
  element_add(QR, Q, R);
  sn_miller(out, pairing->r, in1, curve_x_coord(QR), curve_y_coord(QR));
  sn_miller(e0, pairing->r, in1, curve_x_coord(R), curve_y_coord(R));
  element_invert(e0, e0);
  element_mul(out, out, e0);
  //element_pow_mpz(out, out, p->tateexp);
  element_clear(R);
  element_clear(QR);
}

void pairing_init_singular_with_node(pairing_t pairing, mpz_t q) {
  sn_pairing_data_ptr p;

  mpz_init(pairing->r);
  mpz_sub_ui(pairing->r, q, 1);
  field_init_fp(pairing->Zr, pairing->r);
  pairing->map = sn_pairing;

  p = pairing->data = pbc_malloc(sizeof(sn_pairing_data_t));
  field_init_fp(p->Fq, q);
  field_init_curve_singular_with_node(p->Eq, p->Fq);

  //mpz_init(p->tateexp);
  //mpz_sub_ui(p->tateexp, p->Fq->order, 1);
  //mpz_divexact(p->tateexp, p->tateexp, pairing->r);

  pairing->G2 = pairing->G1 = p->Eq;

  pairing_GT_init(pairing, p->Fq);
}
