#include <stdarg.h>
#include <stdio.h>
#include <stdint.h> // for intptr_t
#include <stdlib.h> //for rand, pbc_malloc, pbc_free
#include <string.h> //for strcmp
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_field.h"
#include "pbc_fp.h"
#include "pbc_param.h"
#include "pbc_pairing.h"
#include "pbc_curve.h"
#include "pbc_random.h"
#include "pbc_memory.h"
#include "pbc_e_param.h"
#include "ecc/param.h"

struct e_param_s {
  mpz_t q;    // Curve is defined over F_q.
  mpz_t r;    // q = h r^2 + 1, r is prime.
  mpz_t h;    // h is 28 times some square.
  mpz_t a, b; // Curve equation is Y^2 = X^3 + aX + b.
  int exp2;
  int exp1;
  int sign1;
  int sign0;
};
typedef struct e_param_s e_param_t[1];
typedef struct e_param_s *e_param_ptr;

struct e_pairing_data_s {
  field_t Fq, Eq;
  int exp2, exp1;
  int sign1, sign0;
  element_t R;
};
typedef struct e_pairing_data_s e_pairing_data_t[1];
typedef struct e_pairing_data_s *e_pairing_data_ptr;

static void e_clear(void *data) {
  e_param_ptr ep = data;
  mpz_clear(ep->q);
  mpz_clear(ep->r);
  mpz_clear(ep->h);
  mpz_clear(ep->a);
  mpz_clear(ep->b);
  pbc_free(data);
}

static void e_out_str(FILE *stream, void *data) {
  e_param_ptr p = data;
  param_out_type(stream, "e");
  param_out_mpz(stream, "q", p->q);
  param_out_mpz(stream, "r", p->r);
  param_out_mpz(stream, "h", p->h);
  param_out_mpz(stream, "a", p->a);
  param_out_mpz(stream, "b", p->b);
  param_out_int(stream, "exp2", p->exp2);
  param_out_int(stream, "exp1", p->exp1);
  param_out_int(stream, "sign1", p->sign1);
  param_out_int(stream, "sign0", p->sign0);
}

static void e_miller_proj(element_t res, element_t P,
    element_ptr QR, element_ptr R,
    e_pairing_data_ptr p) {
  //collate divisions
  int n;
  element_t v, vd;
  element_t v1, vd1;
  element_t Z, Z1;
  element_t a, b, c;
  const element_ptr cca = curve_a_coeff(P);
  element_t e0, e1;
  const element_ptr e2 = a, e3 = b;
  element_t z, z2;
  int i;
  element_ptr Zx, Zy;
  const element_ptr Px = curve_x_coord(P);
  const element_ptr numx = curve_x_coord(QR);
  const element_ptr numy = curve_y_coord(QR);
  const element_ptr denomx = curve_x_coord(R);
  const element_ptr denomy = curve_y_coord(R);

  //convert Z from weighted projective (Jacobian) to affine
  //i.e. (X, Y, Z) --> (X/Z^2, Y/Z^3)
  //also sets z to 1
  #define to_affine() {      \
    element_invert(z, z);    \
    element_square(e0, z);   \
    element_mul(Zx, Zx, e0); \
    element_mul(e0, e0, z);  \
    element_mul(Zy, Zy, e0); \
    element_set1(z);         \
    element_set1(z2);        \
  }

  #define proj_double() {             \
    const element_ptr x = Zx;         \
    const element_ptr y = Zy;         \
    /* e0 = 3x^2 + (cc->a) z^4 */     \
    element_square(e0, x);            \
    /* element_mul_si(e0, e0, 3); */  \
    element_double(e1, e0);           \
    element_add(e0, e0, e1);          \
    element_square(e1, z2);           \
    element_mul(e1, e1, cca);         \
    element_add(e0, e0, e1);          \
                                      \
    /* z_out = 2 y z */               \
    element_mul(z, y, z);             \
    /* element_mul_si(z, z, 2); */    \
    element_double(z, z);             \
    element_square(z2, z);            \
                                      \
    /* e1 = 4 x y^2 */                \
    element_square(e2, y);            \
    element_mul(e1, x, e2);           \
    /* element_mul_si(e1, e1, 4); */  \
    element_double(e1, e1);           \
    element_double(e1, e1);           \
                                      \
    /* x_out = e0^2 - 2 e1 */         \
    /* element_mul_si(e3, e1, 2); */  \
    element_double(e3, e1);           \
    element_square(x, e0);            \
    element_sub(x, x, e3);            \
                                      \
    /* e2 = 8y^4 */                   \
    element_square(e2, e2);           \
    /* element_mul_si(e2, e2, 8); */  \
    element_double(e2, e2);           \
    element_double(e2, e2);           \
    element_double(e2, e2);           \
                                      \
    /* y_out = e0(e1 - x_out) - e2 */ \
    element_sub(e1, e1, x);           \
    element_mul(e0, e0, e1);          \
    element_sub(y, e0, e2);           \
  }

  #define do_tangent(e, edenom) {    \
    /* a = -(3x^2 + cca z^4) */      \
    /* b = 2 y z^3 */                \
    /* c = -(2 y^2 + x a) */         \
    /* a = z^2 a */                  \
    element_square(a, z2);           \
    element_mul(a, a, cca);          \
    element_square(b, Zx);           \
    /* element_mul_si(b, b, 3); */   \
    element_double(e0, b);           \
    element_add(b, b, e0);           \
    element_add(a, a, b);            \
    element_neg(a, a);               \
                                     \
    /* element_mul_si(e0, Zy, 2); */ \
    element_double(e0, Zy);          \
    element_mul(b, e0, z2);          \
    element_mul(b, b, z);            \
                                     \
    element_mul(c, Zx, a);           \
    element_mul(a, a, z2);           \
    element_mul(e0, e0, Zy);         \
    element_add(c, c, e0);           \
    element_neg(c, c);               \
                                     \
    element_mul(e0, a, numx);        \
    element_mul(e1, b, numy);        \
    element_add(e0, e0, e1);         \
    element_add(e0, e0, c);          \
    element_mul(e, e, e0);           \
                                     \
    element_mul(e0, a, denomx);      \
    element_mul(e1, b, denomy);      \
    element_add(e0, e0, e1);         \
    element_add(e0, e0, c);          \
    element_mul(edenom, edenom, e0); \
  }

  #define do_vertical(e, edenom, Ax) { \
    element_mul(e0, numx, z2);         \
    element_sub(e0, e0, Ax);           \
    element_mul(e, e, e0);             \
                                       \
    element_mul(e0, denomx, z2);       \
    element_sub(e0, e0, Ax);           \
    element_mul(edenom, edenom, e0);   \
  }

  #define do_line(e, edenom, A, B) {   \
    element_ptr Ax = curve_x_coord(A); \
    element_ptr Ay = curve_y_coord(A); \
    element_ptr Bx = curve_x_coord(B); \
    element_ptr By = curve_y_coord(B); \
                                       \
    element_sub(b, Bx, Ax);            \
    element_sub(a, Ay, By);            \
    element_mul(c, Ax, By);            \
    element_mul(e0, Ay, Bx);           \
    element_sub(c, c, e0);             \
                                       \
    element_mul(e0, a, numx);          \
    element_mul(e1, b, numy);          \
    element_add(e0, e0, e1);           \
    element_add(e0, e0, c);            \
    element_mul(e, e, e0);             \
                                       \
    element_mul(e0, a, denomx);        \
    element_mul(e1, b, denomy);        \
    element_add(e0, e0, e1);           \
    element_add(e0, e0, c);            \
    element_mul(edenom, edenom, e0);   \
  }

  element_init(a, res->field);
  element_init(b, res->field);
  element_init(c, res->field);
  element_init(e0, res->field);
  element_init(e1, res->field);
  element_init(z, res->field);
  element_init(z2, res->field);
  element_set1(z);
  element_set1(z2);

  element_init(v, res->field);
  element_init(vd, res->field);
  element_init(v1, res->field);
  element_init(vd1, res->field);
  element_init(Z, P->field);
  element_init(Z1, P->field);

  element_set(Z, P);
  Zx = curve_x_coord(Z);
  Zy = curve_y_coord(Z);

  element_set1(v);
  element_set1(vd);
  element_set1(v1);
  element_set1(vd1);

  n = p->exp1;
  for (i=0; i<n; i++) {
    element_square(v, v);
    element_square(vd, vd);
    do_tangent(v, vd);
    proj_double();
    do_vertical(vd, v, Zx);
  }
  to_affine();
  if (p->sign1 < 0) {
    element_set(v1, vd);
    element_set(vd1, v);
    do_vertical(vd1, v1, Zx);
    element_neg(Z1, Z);
  } else {
    element_set(v1, v);
    element_set(vd1, vd);
    element_set(Z1, Z);
  }
  n = p->exp2;
  for (; i<n; i++) {
    element_square(v, v);
    element_square(vd, vd);
    do_tangent(v, vd);
    proj_double();
    do_vertical(vd, v, Zx);
  }
  to_affine();
  element_mul(v, v, v1);
  element_mul(vd, vd, vd1);
  do_line(v, vd, Z, Z1);
  element_add(Z, Z, Z1);
  do_vertical(vd, v, Zx);

  if (p->sign0 > 0) {
    do_vertical(v, vd, Px);
  }

  element_invert(vd, vd);
  element_mul(res, v, vd);

  element_clear(v);
  element_clear(vd);
  element_clear(v1);
  element_clear(vd1);
  element_clear(z);
  element_clear(z2);
  element_clear(Z);
  element_clear(Z1);
  element_clear(a);
  element_clear(b);
  element_clear(c);
  element_clear(e0);
  element_clear(e1);
  #undef to_affine
  #undef proj_double
  #undef do_tangent
  #undef do_vertical
  #undef do_line
}

static void e_miller_affine(element_t res, element_t P,
    element_ptr QR, element_ptr R,
    e_pairing_data_ptr p) {
  //collate divisions
  int n;
  element_t v, vd;
  element_t v1, vd1;
  element_t Z, Z1;
  element_t a, b, c;
  element_t e0, e1;
  const element_ptr Px = curve_x_coord(P);
  const element_ptr cca = curve_a_coeff(P);
  element_ptr Zx, Zy;
  int i;
  const element_ptr numx = curve_x_coord(QR);
  const element_ptr numy = curve_y_coord(QR);
  const element_ptr denomx = curve_x_coord(R);
  const element_ptr denomy = curve_y_coord(R);

  #define do_vertical(e, edenom, Ax) { \
    element_sub(e0, numx, Ax);         \
    element_mul(e, e, e0);             \
                                       \
    element_sub(e0, denomx, Ax);       \
    element_mul(edenom, edenom, e0);   \
  }

  #define do_tangent(e, edenom) {                      \
    /* a = -slope_tangent(A.x, A.y); */                \
    /* b = 1; */                                       \
    /* c = -(A.y + a * A.x); */                        \
    /* but we multiply by 2*A.y to avoid division */   \
                                                       \
    /* a = -Ax * (Ax + Ax + Ax + twicea_2) - a_4; */   \
    /* Common curves: a2 = 0 (and cc->a is a_4), so */ \
    /* a = -(3 Ax^2 + cc->a) */                        \
    /* b = 2 * Ay */                                   \
    /* c = -(2 Ay^2 + a Ax); */                        \
                                                       \
    element_square(a, Zx);                             \
    element_mul_si(a, a, 3);                           \
    element_add(a, a, cca);                            \
    element_neg(a, a);                                 \
                                                       \
    element_add(b, Zy, Zy);                            \
                                                       \
    element_mul(e0, b, Zy);                            \
    element_mul(c, a, Zx);                             \
    element_add(c, c, e0);                             \
    element_neg(c, c);                                 \
                                                       \
    element_mul(e0, a, numx);                          \
    element_mul(e1, b, numy);                          \
    element_add(e0, e0, e1);                           \
    element_add(e0, e0, c);                            \
    element_mul(e, e, e0);                             \
                                                       \
    element_mul(e0, a, denomx);                        \
    element_mul(e1, b, denomy);                        \
    element_add(e0, e0, e1);                           \
    element_add(e0, e0, c);                            \
    element_mul(edenom, edenom, e0);                   \
  }

  #define do_line(e, edenom, A, B) {   \
    element_ptr Ax = curve_x_coord(A); \
    element_ptr Ay = curve_y_coord(A); \
    element_ptr Bx = curve_x_coord(B); \
    element_ptr By = curve_y_coord(B); \
                                       \
    element_sub(b, Bx, Ax);            \
    element_sub(a, Ay, By);            \
    element_mul(c, Ax, By);            \
    element_mul(e0, Ay, Bx);           \
    element_sub(c, c, e0);             \
                                       \
    element_mul(e0, a, numx);          \
    element_mul(e1, b, numy);          \
    element_add(e0, e0, e1);           \
    element_add(e0, e0, c);            \
    element_mul(e, e, e0);             \
                                       \
    element_mul(e0, a, denomx);        \
    element_mul(e1, b, denomy);        \
    element_add(e0, e0, e1);           \
    element_add(e0, e0, c);            \
    element_mul(edenom, edenom, e0);   \
  }

  element_init(a, res->field);
  element_init(b, res->field);
  element_init(c, res->field);
  element_init(e0, res->field);
  element_init(e1, res->field);

  element_init(v, res->field);
  element_init(vd, res->field);
  element_init(v1, res->field);
  element_init(vd1, res->field);
  element_init(Z, P->field);
  element_init(Z1, P->field);

  element_set(Z, P);
  Zx = curve_x_coord(Z);
  Zy = curve_y_coord(Z);

  element_set1(v);
  element_set1(vd);
  element_set1(v1);
  element_set1(vd1);

  n = p->exp1;
  for (i=0; i<n; i++) {
    element_square(v, v);
    element_square(vd, vd);
    do_tangent(v, vd);
    element_double(Z, Z);
    do_vertical(vd, v, Zx);
  }
  if (p->sign1 < 0) {
    element_set(v1, vd);
    element_set(vd1, v);
    do_vertical(vd1, v1, Zx);
    element_neg(Z1, Z);
  } else {
    element_set(v1, v);
    element_set(vd1, vd);
    element_set(Z1, Z);
  }
  n = p->exp2;
  for (; i<n; i++) {
    element_square(v, v);
    element_square(vd, vd);
    do_tangent(v, vd);
    element_double(Z, Z);
    do_vertical(vd, v, Zx);
  }
  element_mul(v, v, v1);
  element_mul(vd, vd, vd1);
  do_line(v, vd, Z, Z1);
  element_add(Z, Z, Z1);
  do_vertical(vd, v, Zx);

  if (p->sign0 > 0) {
    do_vertical(v, vd, Px);
  }

  element_invert(vd, vd);
  element_mul(res, v, vd);

  element_clear(v);
  element_clear(vd);
  element_clear(v1);
  element_clear(vd1);
  element_clear(Z);
  element_clear(Z1);
  element_clear(a);
  element_clear(b);
  element_clear(c);
  element_clear(e0);
  element_clear(e1);
  #undef do_vertical
  #undef do_tangent
  #undef do_line
}

static void (*e_miller_fn)(element_t res, element_t P,
    element_ptr QR, element_ptr R,
    e_pairing_data_ptr p);

static void e_pairing(element_ptr out, element_ptr in1, element_ptr in2,
    pairing_t pairing) {
  e_pairing_data_ptr p = pairing->data;
  element_ptr Q = in2;
  element_t QR;
  element_init(QR, p->Eq);
  element_add(QR, Q, p->R);
  e_miller_fn(out, in1, QR, p->R, p);
  element_pow_mpz(out, out, pairing->phikonr);
  element_clear(QR);
}

// in1, in2 are from E(F_q), out from F_q^2.
// Pairing via elliptic nets (see Stange).
static void e_pairing_ellnet(element_ptr out, element_ptr in1, element_ptr in2,
    pairing_t pairing) {
  const element_ptr a = curve_a_coeff(in1);
  const element_ptr b = curve_b_coeff(in1);

  element_ptr x = curve_x_coord(in1);
  element_ptr y = curve_y_coord(in1);

  element_ptr x2 = curve_x_coord(in2);
  element_ptr y2 = curve_y_coord(in2);

  //notation: cmi means c_{k-i}, ci means c_{k+i}
  element_t cm3, cm2, cm1, c0, c1, c2, c3, c4;
  element_t dm1, d0, d1;
  element_t A, B, C;

  element_init_same_as(cm3, x);
  element_init_same_as(cm2, x);
  element_init_same_as(cm1, x);
  element_init_same_as(c0, x);
  element_init_same_as(c1, x);
  element_init_same_as(c2, x);
  element_init_same_as(c3, x);
  element_init_same_as(c4, x);
  element_init_same_as(C, x);

  element_init_same_as(dm1, out);
  element_init_same_as(d0, out);
  element_init_same_as(d1, out);
  element_init_same_as(A, x);
  element_init_same_as(B, out);

  // c1 = 2y
  // cm3 = -2y
  element_double(c1, y);
  element_neg(cm3, c1);

  //use c0, cm1, cm2, C, c4 as temp variables for now
  //compute c3, c2
  element_square(cm2, x);
  element_square(C, cm2);
  element_mul(cm1, b, x);
  element_double(cm1, cm1);
  element_square(c4, a);

  element_mul(c2, cm1, cm2);
  element_double(c2, c2);
  element_mul(c0, a, C);
  element_add(c2, c2, c0);
  element_mul(c0, c4, cm2);
  element_sub(c2, c2, c0);
  element_double(c0, c2);
  element_double(c0, c0);
  element_add(c2, c2, c0);

  element_mul(c0, cm1, a);
  element_square(c3, b);
  element_double(c3, c3);
  element_double(c3, c3);
  element_add(c0, c0, c3);
  element_double(c0, c0);
  element_mul(c3, a, c4);
  element_add(c0, c0, c3);
  element_sub(c2, c2, c0);
  element_mul(c0, cm2, C);
  element_add(c3, c0, c2);
  element_mul(c3, c3, c1);
  element_double(c3, c3);

  element_mul(c0, a, cm2);
  element_add(c0, c0, cm1);
  element_double(c0, c0);
  element_add(c0, c0, C);
  element_double(c2, c0);
  element_add(c0, c0, c2);
  element_sub(c2, c0, c4);

  // c0 = 1
  // cm2 = -1
  element_set1(c0);
  element_neg(cm2, c0);

  // c4 = c_5 = c_2^3 c_4 - c_3^3 = c1^3 c3 - c2^3
  element_square(C, c1);
  element_mul(c4, C, c1);
  element_mul(c4, c4, c3);
  element_square(C, c2);
  element_mul(C, C, c2);
  element_sub(c4, c4, C);

  //compute A, B, d1 (which is d_2 since k = 1)
  element_sub(A, x, x2);
  element_double(C, x);
  element_add(C, C, x2);
  element_square(cm1, A);
  element_mul(cm1, C, cm1);
  element_add(d1, y, y2);
  element_square(d1, d1);
  element_sub(B, cm1, d1);
  element_invert(B, B);
  element_invert(A, A);

  element_sub(d1, y, y2);
  element_mul(d1, d1, A);
  element_square(d1, d1);
  element_sub(d1, C, d1);

  // cm1 = 0
  // C = (2y)^-1
  element_set0(cm1);
  element_invert(C, c1);

  element_set1(dm1);
  element_set1(d0);

  element_t sm2, sm1;
  element_t s0, s1, s2, s3;
  element_t tm2, tm1;
  element_t t0, t1, t2, t3;
  element_t e0, e1;
  element_t u, v;

  element_init_same_as(sm2, x);
  element_init_same_as(sm1, x);
  element_init_same_as(s0, x);
  element_init_same_as(s1, x);
  element_init_same_as(s2, x);
  element_init_same_as(s3, x);

  element_init_same_as(tm2, x);
  element_init_same_as(tm1, x);
  element_init_same_as(t0, x);
  element_init_same_as(t1, x);
  element_init_same_as(t2, x);
  element_init_same_as(t3, x);

  element_init_same_as(e0, x);
  element_init_same_as(e1, x);

  element_init_same_as(u, d0);
  element_init_same_as(v, d0);

  int m = mpz_sizeinbase(pairing->r, 2) - 2;
  for (;;) {
    element_square(sm2, cm2);
    element_square(sm1, cm1);
    element_square(s0, c0);
    element_square(s1, c1);
    element_square(s2, c2);
    element_square(s3, c3);

    element_mul(tm2, cm3, cm1);
    element_mul(tm1, cm2, c0);
    element_mul(t0, cm1, c1);
    element_mul(t1, c0, c2);
    element_mul(t2, c1, c3);
    element_mul(t3, c2, c4);

    element_square(u, d0);
    element_mul(v, dm1, d1);

    if (mpz_tstbit(pairing->r, m)) {
      //double-and-add
      element_mul(e0, t0, sm2);
      element_mul(e1, tm2, s0);
      element_sub(cm3, e0, e1);
      element_mul(cm3, cm3, C);

      element_mul(e0, t0, sm1);
      element_mul(e1, tm1, s0);
      element_sub(cm2, e0, e1);

      element_mul(e0, t1, sm1);
      element_mul(e1, tm1, s1);
      element_sub(cm1, e0, e1);
      element_mul(cm1, cm1, C);

      element_mul(e0, t1, s0);
      element_mul(e1, t0, s1);
      element_sub(c0, e0, e1);

      element_mul(e0, t2, s0);
      element_mul(e1, t0, s2);
      element_sub(c1, e0, e1);
      element_mul(c1, c1, C);

      element_mul(e0, t2, s1);
      element_mul(e1, t1, s2);
      element_sub(c2, e0, e1);

      element_mul(e0, t3, s1);
      element_mul(e1, t1, s3);
      element_sub(c3, e0, e1);
      element_mul(c3, c3, C);

      element_mul(e0, t3, s2);
      element_mul(e1, t2, s3);
      element_sub(c4, e0, e1);

      element_mul(out, u, t0);
      element_mul(dm1, v, s0);
      element_sub(dm1, dm1, out);

      element_mul(out, u, t1);
      element_mul(d0, v, s1);
      element_sub(d0, d0, out);
      element_mul(d0, d0, A);

      element_mul(out, u, t2);
      element_mul(d1, v, s2);
      element_sub(d1, d1, out);
      element_mul(d1, d1, B);
    } else {
      //double
      element_mul(e0, tm1, sm2);
      element_mul(e1, tm2, sm1);
      element_sub(cm3, e0, e1);

      element_mul(e0, t0, sm2);
      element_mul(e1, tm2, s0);
      element_sub(cm2, e0, e1);
      element_mul(cm2, cm2, C);

      element_mul(e0, t0, sm1);
      element_mul(e1, tm1, s0);
      element_sub(cm1, e0, e1);

      element_mul(e0, t1, sm1);
      element_mul(e1, tm1, s1);
      element_sub(c0, e0, e1);
      element_mul(c0, c0, C);

      element_mul(e0, t1, s0);
      element_mul(e1, t0, s1);
      element_sub(c1, e0, e1);

      element_mul(e0, t2, s0);
      element_mul(e1, t0, s2);
      element_sub(c2, e0, e1);
      element_mul(c2, c2, C);

      element_mul(e0, t2, s1);
      element_mul(e1, t1, s2);
      element_sub(c3, e0, e1);

      element_mul(e0, t3, s1);
      element_mul(e1, t1, s3);
      element_sub(c4, e0, e1);
      element_mul(c4, c4, C);

      element_mul(out, u, tm1);
      element_mul(dm1, v, sm1);
      element_sub(dm1, dm1, out);

      element_mul(out, u, t0);
      element_mul(d0, v, s0);
      element_sub(d0, d0, out);

      element_mul(out, u, t1);
      element_mul(d1, v, s1);
      element_sub(d1, d1, out);
      element_mul(d1, d1, A);
    }
    if (!m) break;
    m--;
  }
  element_invert(c1, c1);
  element_mul(d1, d1, c1);

  element_pow_mpz(out, d1, pairing->phikonr);

  element_clear(dm1);
  element_clear(d0);
  element_clear(d1);

  element_clear(cm3);
  element_clear(cm2);
  element_clear(cm1);
  element_clear(c0);
  element_clear(c1);
  element_clear(c2);
  element_clear(c3);
  element_clear(c4);

  element_clear(sm2);
  element_clear(sm1);
  element_clear(s0);
  element_clear(s1);
  element_clear(s2);
  element_clear(s3);

  element_clear(tm2);
  element_clear(tm1);
  element_clear(t0);
  element_clear(t1);
  element_clear(t2);
  element_clear(t3);

  element_clear(e0);
  element_clear(e1);
  element_clear(A);
  element_clear(B);
  element_clear(C);
  element_clear(u);
  element_clear(v);
}

static void phi_identity(element_ptr out, element_ptr in, pairing_ptr pairing) {
  (void) pairing;
  element_set(out, in);
}

static void e_pairing_option_set(pairing_t pairing, char *key, char *value) {
  //TODO: this affects every type E pairing!
  UNUSED_VAR(pairing);
  if (!strcmp(key, "method")) {
    if (!strcmp(value, "miller")) {
      pairing->map = e_pairing;
      e_miller_fn = e_miller_proj;
    } else if (!strcmp(value, "miller-affine")) {
      pairing->map = e_pairing;
      e_miller_fn = e_miller_affine;
    } else if (!strcmp(value, "shipsey-stange")) {
      pairing->map = e_pairing_ellnet;
    }
  }
}

static void e_pairing_clear(pairing_t pairing) {
  field_clear(pairing->GT);
  e_pairing_data_ptr p = pairing->data;
  field_clear(p->Fq);
  field_clear(p->Eq);
  element_clear(p->R);
  pbc_free(p);

  mpz_clear(pairing->phikonr);
  mpz_clear(pairing->r);
  field_clear(pairing->Zr);
}

static void e_finalpow(element_ptr e) {
  element_pow_mpz(e->data, e->data, e->field->pairing->phikonr);
}

static void e_init_pairing(pairing_t pairing, void *data) {
  e_param_ptr param = data;
  e_pairing_data_ptr p;
  element_t a, b;

  mpz_init(pairing->r);
  mpz_set(pairing->r, param->r);
  field_init_fp(pairing->Zr, pairing->r);
  pairing->map = e_pairing;
  e_miller_fn = e_miller_proj;

  p = pairing->data = pbc_malloc(sizeof(e_pairing_data_t));
  p->exp2 = param->exp2;
  p->exp1 = param->exp1;
  p->sign1 = param->sign1;
  p->sign0 = param->sign0;
  field_init_fp(p->Fq, param->q);
  element_init(a, p->Fq);
  element_init(b, p->Fq);
  element_set_mpz(a, param->a);
  element_set_mpz(b, param->b);
  field_init_curve_ab(p->Eq, a, b, pairing->r, param->h);

  //k=1, hence phikonr = (p-1)/r
  mpz_init(pairing->phikonr);
  mpz_sub_ui(pairing->phikonr, p->Fq->order, 1);
  mpz_divexact(pairing->phikonr, pairing->phikonr, pairing->r);

  pairing->G2 = pairing->G1 = p->Eq;
  pairing_GT_init(pairing, p->Fq);
  pairing->finalpow = e_finalpow;
  pairing->phi = phi_identity;
  pairing->option_set = e_pairing_option_set;
  pairing->clear_func = e_pairing_clear;

  element_init(p->R, p->Eq);
  curve_set_gen_no_cofac(p->R);

  element_clear(a);
  element_clear(b);
}

static void e_init(pbc_param_ptr p) {
  static pbc_param_interface_t interface = {{
    e_clear,
    e_init_pairing,
    e_out_str,
  }};
  p->api = interface;
  e_param_ptr ep = p->data = pbc_malloc(sizeof(*ep));
  mpz_init(ep->q);
  mpz_init(ep->r);
  mpz_init(ep->h);
  mpz_init(ep->a);
  mpz_init(ep->b);
}

// Public interface:

int pbc_param_init_e(pbc_param_ptr par, struct symtab_s *tab) {
  e_init(par);
  e_param_ptr p = par->data;

  int err = 0;
  err += lookup_mpz(p->q, tab, "q");
  err += lookup_mpz(p->r, tab, "r");
  err += lookup_mpz(p->h, tab, "h");
  err += lookup_mpz(p->a, tab, "a");
  err += lookup_mpz(p->b, tab, "b");
  err += lookup_int(&p->exp2, tab, "exp2");
  err += lookup_int(&p->exp1, tab, "exp1");
  err += lookup_int(&p->sign1, tab, "sign1");
  err += lookup_int(&p->sign0, tab, "sign0");
  return err;
}

void pbc_param_init_e_gen(pbc_param_t par, int rbits, int qbits) {
  e_init(par);
  e_param_ptr p = par->data;
  //3 takes 2 bits to represent
  int hbits = (qbits - 2) / 2 - rbits;
  mpz_ptr q = p->q;
  mpz_ptr r = p->r;
  mpz_ptr h = p->h;
  mpz_t n;
  field_t Fq;
  field_t cc;
  element_t j;
  int found = 0;

  //won't find any curves is hbits is too low
  if (hbits < 3) hbits = 3;

  mpz_init(n);

  do {
    int i;
    mpz_set_ui(r, 0);

    if (rand() % 2) {
      p->exp2 = rbits - 1;
      p->sign1 = 1;
    } else {
      p->exp2 = rbits;
      p->sign1 = -1;
    }
    mpz_setbit(r, p->exp2);

    p->exp1 = (rand() % (p->exp2 - 1)) + 1;
    //use q as a temp variable
    mpz_set_ui(q, 0);
    mpz_setbit(q, p->exp1);

    if (p->sign1 > 0) {
      mpz_add(r, r, q);
    } else {
      mpz_sub(r, r, q);
    }

    if (rand() % 2) {
      p->sign0 = 1;
      mpz_add_ui(r, r, 1);
    } else {
      p->sign0 = -1;
      mpz_sub_ui(r, r, 1);
    }
    if (!mpz_probab_prime_p(r, 10)) continue;
    for (i=0; i<10; i++) {
      //use q as a temp variable
      mpz_set_ui(q, 0);
      mpz_setbit(q, hbits + 1);
      pbc_mpz_random(h, q);
      mpz_mul(h, h, h);
      mpz_mul_ui(h, h, 3);
      //finally q takes the value it should
      mpz_mul(n, r, r);
      mpz_mul(n, n, h);
      mpz_add_ui(q, n, 1);
      if (mpz_probab_prime_p(q, 10)) {
        found = 1;
        break;
      }
    }
  } while (!found);
  /*
  do {
    mpz_set_ui(r, 0);
    mpz_setbit(r, rbits);
    pbc_mpz_random(r, r);
    mpz_nextprime(r, r);
    mpz_mul(n, r, r);
    mpz_mul_ui(n, n, 3);
    mpz_add_ui(q, n, 1);
  } while (!mpz_probab_prime_p(q, 10));
  */

  field_init_fp(Fq, q);
  element_init(j, Fq);
  element_set_si(j, 1);
  field_init_curve_b(cc, j, n, NULL);
  element_clear(j);
  // We may need to twist it.
  {
    // Pick a random point P and twist the curve if P has the wrong order.
    element_t P;
    element_init(P, cc);
    element_random(P);
    element_mul_mpz(P, P, n);
    if (!element_is0(P)) field_reinit_curve_twist(cc);
    element_clear(P);
  }
  element_to_mpz(p->a, curve_field_a_coeff(cc));
  element_to_mpz(p->b, curve_field_b_coeff(cc));

  mpz_clear(n);
}
