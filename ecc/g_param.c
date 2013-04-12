#include <stdarg.h>
#include <stdio.h>
#include <stdint.h> // for intptr_t
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_field.h"
#include "pbc_poly.h"
#include "pbc_hilbert.h"
#include "pbc_fp.h"
#include "pbc_fieldquadratic.h"
#include "pbc_mnt.h"
#include "pbc_curve.h"
#include "pbc_param.h"
#include "pbc_pairing.h"
#include "pbc_memory.h"
#include "pbc_g_param.h"
#include "ecc/param.h"

struct g_param_s {
  mpz_t q;    // Curve defined over F_q.
  mpz_t n;    // n = #E(F_q) (= q - t + 1)
  mpz_t h;    // h * r = n, r is prime
  mpz_t r;
  mpz_t a, b; // E: y^2 = x^3 + ax + b

  // k = 10 for these curves.
  mpz_t nk;      // #E(F_q^k)
  mpz_t hk;      // hk * r^2 = nk
  mpz_t *coeff;  //Coefficients of polynomial used to extend F_q by k/2
  mpz_t nqr;     // Quadratic nonresidue in F_q^d that lies in F_q.
};

typedef struct g_param_s g_param_t[1];
typedef struct g_param_s *g_param_ptr;

struct mnt_pairing_data_s {
  field_t Fq, Fqx, Fqd, Fqk;
  field_t Eq, Etwist;
  element_t nqrinv, nqrinv2;
  element_t xpowq, xpowq2, xpowq3, xpowq4;
};
typedef struct mnt_pairing_data_s mnt_pairing_data_t[1];
typedef struct mnt_pairing_data_s *mnt_pairing_data_ptr;

static void g_clear(void *data) {
  g_param_ptr param = data;
  int i;
  mpz_clear(param->q);
  mpz_clear(param->n);
  mpz_clear(param->h);
  mpz_clear(param->r);
  mpz_clear(param->a);
  mpz_clear(param->b);
  mpz_clear(param->nk);
  mpz_clear(param->hk);
  mpz_clear(param->nqr);
  for (i = 0; i < 5; i++) {
    mpz_clear(param->coeff[i]);
  }
  pbc_free(param->coeff);
  pbc_free(data);
}

static void g_out_str(FILE *stream, void *data) {
  g_param_ptr p = data;
  int i;
  char s[8];
  param_out_type(stream, "g");
  param_out_mpz(stream, "q", p->q);
  param_out_mpz(stream, "n", p->n);
  param_out_mpz(stream, "h", p->h);
  param_out_mpz(stream, "r", p->r);
  param_out_mpz(stream, "a", p->a);
  param_out_mpz(stream, "b", p->b);
  param_out_mpz(stream, "nk", p->nk);
  param_out_mpz(stream, "hk", p->hk);
  for (i=0; i<5; i++) {
    sprintf(s, "coeff%d", i);
    param_out_mpz(stream, s, p->coeff[i]);
  }
  param_out_mpz(stream, "nqr", p->nqr);
}

static inline void d_miller_evalfn(element_t e0,
    element_t a, element_t b, element_t c,
    element_t Qx, element_t Qy) {
  //a, b, c are in Fq
  //point Q is (Qx, Qy * sqrt(nqr)) where nqr is used to construct
  //the quadratic field extension Fqk of Fqd
  element_ptr re_out = element_x(e0);
  element_ptr im_out = element_y(e0);

  int i;
  int d = polymod_field_degree(re_out->field);
  for (i=0; i<d; i++) {
    element_mul(element_item(re_out, i), element_item(Qx, i), a);
    element_mul(element_item(im_out, i), element_item(Qy, i), b);
  }
  element_add(element_item(re_out, 0), element_item(re_out, 0), c);
}

static void cc_miller_no_denom_proj(element_t res, mpz_t q, element_t P,
    element_ptr Qx, element_ptr Qy) {
  int m;
  element_t v;
  element_t Z;
  element_t a, b, c;
  element_t t0, t1;
  element_ptr t2 = a, t3 = b, t4 = c;
  element_t e0;
  element_t z, z2;
  element_ptr Zx, Zy;
  const element_ptr curve_a = curve_a_coeff(P);
  const element_ptr Px = curve_x_coord(P);
  const element_ptr Py = curve_y_coord(P);

  #define proj_double() {             \
    /* t0 = 3x^2 + (curve_a) z^4 */   \
    element_square(t0, Zx);           \
    /* element_mul_si(t0, t0, 3); */  \
    element_double(t1, t0);           \
    element_add(t0, t0, t1);          \
    element_square(t1, z2);           \
    element_mul(t1, t1, curve_a);     \
    element_add(t0, t0, t1);          \
                                      \
    /* z_out = 2 y z */               \
    element_mul(z, Zy, z);            \
    /* element_mul_si(z, z, 2); */    \
    element_double(z, z);             \
    element_square(z2, z);            \
                                      \
    /* t1 = 4 x y^2 */                \
    element_square(t2, Zy);           \
    element_mul(t1, Zx, t2);          \
    /* element_mul_si(t1, t1, 4); */  \
    element_double(t1, t1);           \
    element_double(t1, t1);           \
                                      \
    /* x_out = t0^2 - 2 t1 */         \
    /* element_mul_si(t3, t1, 2); */  \
    element_double(t3, t1);           \
    element_square(Zx, t0);           \
    element_sub(Zx, Zx, t3);          \
                                      \
    /* t2 = 8y^4 */                   \
    element_square(t2, t2);           \
    /* element_mul_si(t2, t2, 8); */  \
    element_double(t2, t2);           \
    element_double(t2, t2);           \
    element_double(t2, t2);           \
                                      \
    /* y_out = t0(t1 - x_out) - t2 */ \
    element_sub(t1, t1, Zx);          \
    element_mul(t0, t0, t1);          \
    element_sub(Zy, t0, t2);          \
  }

  #define proj_mixin() {                        \
    /* t2 = Px z^2 */                           \
    element_mul(t2, z2, Px);                    \
                                                \
    /* t3 = Zx - t2 */                          \
    element_sub(t3, Zx, t2);                    \
                                                \
    /* t0 = Py z^3 */                           \
    element_mul(t0, z2, Py);                    \
    element_mul(t0, t0, z);                     \
                                                \
    /* t1 = Zy - t0 */                          \
    element_sub(t1, Zy, t0);                    \
                                                \
    /* e7 = Zx + t2, use t2 to double for e7 */ \
    element_add(t2, Zx, t2);                    \
                                                \
    /* e8 = Zy + t0, use t0 to double for e8 */ \
    element_add(t0, Zy, t0);                    \
                                                \
    /* z = z t3 */                              \
    element_mul(z, z, t3);                      \
    element_square(z2, z);                      \
                                                \
    /* Zx = t1^2 - e7 t3^2 */                   \
    /* t3 now holds t3^3, */                    \
    /* t4 holds e7 t3^2 */                      \
    element_square(t4, t3);                     \
    element_mul(t3, t4, t3);                    \
    element_square(Zx, t1);                     \
    element_mul(t4, t2, t4);                    \
    element_sub(Zx, Zx, t4);                    \
                                                \
    /* t4 = e7 t3^2 - 2 Zx */                   \
    element_sub(t4, t4, Zx);                    \
    element_sub(t4, t4, Zx);                    \
                                                \
    /* Zy = (t4 t1 - e8 t3^3)/2 */              \
    element_mul(t4, t4, t1);                    \
    element_mul(t0, t0, t3);                    \
    element_sub(t4, t4, t0);                    \
    element_halve(Zy, t4);                      \
  }

  #define do_tangent() {                  \
    /* a = -(3x^2 + cca z^4) */           \
    /* b = 2 y z^3 */                     \
    /* c = -(2 y^2 + x a) */              \
    /* a = z^2 a */                       \
    element_square(a, z2);                \
    element_mul(a, a, curve_a);           \
    element_square(b, Zx);                \
    /* element_mul_si(b, b, 3); */        \
    element_double(t0, b);                \
    element_add(b, b, t0);                \
    element_add(a, a, b);                 \
    element_neg(a, a);                    \
                                          \
    element_mul(b, z, z2);                \
    element_mul(b, b, Zy);                \
    element_mul_si(b, b, 2);              \
                                          \
    element_mul(c, Zx, a);                \
    element_mul(a, a, z2);                \
    element_square(t0, Zy);               \
    element_mul_si(t0, t0, 2);            \
    element_add(c, c, t0);                \
    element_neg(c, c);                    \
                                          \
    d_miller_evalfn(e0, a, b, c, Qx, Qy); \
    element_mul(v, v, e0);                \
  }

  #define do_line() {                     \
    /* a = -(Py z^3 - Zy) */              \
    /* b = Px z^3 - Zx z */               \
    /* c = Zx z Py - Zy Px; */            \
                                          \
    element_mul(t0, Zx, z);               \
    element_mul(t1, z2, z);               \
                                          \
    element_mul(a, Py, t1);               \
    element_sub(a, Zy, a);                \
                                          \
    element_mul(b, Px, t1);               \
    element_sub(b, b, t0);                \
                                          \
    element_mul(t0, t0, Py);              \
    element_mul(c, Zy, Px);               \
    element_sub(c, t0, c);                \
                                          \
    d_miller_evalfn(e0, a, b, c, Qx, Qy); \
    element_mul(v, v, e0);                \
  }

  element_init(a, Px->field);
  element_init(b, a->field);
  element_init(c, a->field);
  element_init(t0, a->field);
  element_init(t1, a->field);
  element_init(e0, res->field);
  element_init(z, a->field);
  element_init(z2, a->field);
  element_set1(z);
  element_set1(z2);

  element_init(v, res->field);
  element_init(Z, P->field);

  element_set(Z, P);
  Zx = curve_x_coord(Z);
  Zy = curve_x_coord(Z);

  element_set1(v);
  m = mpz_sizeinbase(q, 2) - 2;

  for(;;) {
    do_tangent();
    if (!m) break;
    proj_double();
    if (mpz_tstbit(q, m)) {
      do_line();
      proj_mixin();
    }
    m--;
    element_square(v, v);
  }

  element_set(res, v);

  element_clear(v);
  element_clear(Z);
  element_clear(a);
  element_clear(b);
  element_clear(c);
  element_clear(t0);
  element_clear(t1);
  element_clear(e0);
  element_clear(z);
  element_clear(z2);
  #undef proj_double
  #undef proj_mixin
  #undef do_tangent
  #undef do_line
}

static void cc_miller_no_denom_affine(element_t res, mpz_t q, element_t P,
    element_ptr Qx, element_ptr Qy) {
  int m;
  element_t v;
  element_t Z;
  element_t a, b, c;
  element_t t0;
  element_t e0;
  const element_ptr cca = curve_a_coeff(P);
  const element_ptr Px = curve_x_coord(P);
  const element_ptr Py = curve_y_coord(P);
  element_ptr Zx, Zy;

  /* TODO: when exactly is this not needed?
  void do_vertical(void)
  {
    mapbase(e0, Z->x);
    element_sub(e0, Qx, e0);
    element_mul(v, v, e0);
  }
  */

  #define do_tangent() {                  \
    /* a = -(3 Zx^2 + cc->a) */           \
    /* b = 2 * Zy */                      \
    /* c = -(2 Zy^2 + a Zx); */           \
    element_square(a, Zx);                \
    element_mul_si(a, a, 3);              \
    element_add(a, a, cca);               \
    element_neg(a, a);                    \
                                          \
    element_add(b, Zy, Zy);               \
                                          \
    element_mul(t0, b, Zy);               \
    element_mul(c, a, Zx);                \
    element_add(c, c, t0);                \
    element_neg(c, c);                    \
                                          \
    d_miller_evalfn(e0, a, b, c, Qx, Qy); \
    element_mul(v, v, e0);                \
  }

  #define do_line() {                     \
    /* a = -(B.y - A.y) / (B.x - A.x); */ \
    /* b = 1; */                          \
    /* c = -(A.y + a * A.x); */           \
    /* but we'll multiply by B.x - A.x */ \
    /* to avoid division */               \
                                          \
    element_sub(b, Px, Zx);               \
    element_sub(a, Zy, Py);               \
    element_mul(t0, b, Zy);               \
    element_mul(c, a, Zx);                \
    element_add(c, c, t0);                \
    element_neg(c, c);                    \
                                          \
    d_miller_evalfn(e0, a, b, c, Qx, Qy); \
    element_mul(v, v, e0);                \
  }

  element_init(a, Px->field);
  element_init(b, a->field);
  element_init(c, a->field);
  element_init(t0, a->field);
  element_init(e0, res->field);

  element_init(v, res->field);
  element_init(Z, P->field);

  element_set(Z, P);
  Zx = curve_x_coord(Z);
  Zy = curve_y_coord(Z);

  element_set1(v);
  m = mpz_sizeinbase(q, 2) - 2;

  for(;;) {
    do_tangent();
    if (!m) break;
    element_double(Z, Z);
    if (mpz_tstbit(q, m)) {
      do_line();
      element_add(Z, Z, P);
    }
    m--;
    element_square(v, v);
  }

  element_set(res, v);

  element_clear(v);
  element_clear(Z);
  element_clear(a);
  element_clear(b);
  element_clear(c);
  element_clear(t0);
  element_clear(e0);
  #undef do_tangent
  #undef do_line
}

// Requires cofactor is even.
// Requires in != out.
// Mangles in.
static void lucas_even(element_ptr out, element_ptr in, mpz_t cofactor) {
  element_t temp;
  element_init_same_as(temp, out);
  element_ptr in0 = element_x(in);
  element_ptr in1 = element_y(in);
  element_ptr v0 = element_x(out);
  element_ptr v1 = element_y(out);
  element_ptr t0 = element_x(temp);
  element_ptr t1 = element_y(temp);
  int j;

  element_set_si(t0, 2);
  element_double(t1, in0);

  element_set(v0, t0);
  element_set(v1, t1);

  j = mpz_sizeinbase(cofactor, 2) - 1;
  for (;;) {
    if (!j) {
      element_mul(v1, v0, v1);
      element_sub(v1, v1, t1);
      element_square(v0, v0);
      element_sub(v0, v0, t0);
      break;
    }
    if (mpz_tstbit(cofactor, j)) {
      element_mul(v0, v0, v1);
      element_sub(v0, v0, t1);
      element_square(v1, v1);
      element_sub(v1, v1, t0);
    } else {
      element_mul(v1, v0, v1);
      element_sub(v1, v1, t1);
      element_square(v0, v0);
      element_sub(v0, v0, t0);
    }
    j--;
  }

  //assume cofactor = (q^2 - q + 1) / r is odd
  //thus v1 = V_k, v0 = V_{k-1}
  //   U = (P v1 - 2 v0) / (P^2 - 4)

  element_double(v0, v0);
  element_mul(in0, t1, v1);
  element_sub(in0, in0, v0);

  element_square(t1, t1);
  element_sub(t1, t1, t0);
  element_sub(t1, t1, t0);

  element_halve(v0, v1);
  element_div(v1, in0, t1);
  element_mul(v1, v1, in1);
  element_clear(temp);
}

static void tatepower10(element_ptr out, element_ptr in, pairing_t pairing) {
  mnt_pairing_data_ptr p = pairing->data;
  element_t e0, e1, e2, e3;
  element_init(e0, p->Fqk);
  element_init(e1, p->Fqd);
  element_init(e2, p->Fqd);
  element_init(e3, p->Fqk);
  element_ptr e0re = element_x(e0);
  element_ptr e0im = element_y(e0);
  element_ptr e0re0 = ((element_t *) e0re->data)[0];
  element_ptr e0im0 = ((element_t *) e0im->data)[0];
  element_t *inre = element_x(in)->data;
  element_t *inim = element_y(in)->data;
  //see thesis
  #define qpower(sign) {                         \
    polymod_const_mul(e2, inre[1], p->xpowq);    \
    element_set(e0re, e2);                       \
    polymod_const_mul(e2, inre[2], p->xpowq2);   \
    element_add(e0re, e0re, e2);                 \
    polymod_const_mul(e2, inre[3], p->xpowq3);   \
    element_add(e0re, e0re, e2);                 \
    polymod_const_mul(e2, inre[4], p->xpowq4);   \
    element_add(e0re, e0re, e2);                 \
    element_add(e0re0, e0re0, inre[0]);          \
                                                 \
    if (sign > 0) {                              \
      polymod_const_mul(e2, inim[1], p->xpowq);  \
      element_set(e0im, e2);                     \
      polymod_const_mul(e2, inim[2], p->xpowq2); \
      element_add(e0im, e0im, e2);               \
      polymod_const_mul(e2, inim[3], p->xpowq3); \
      element_add(e0im, e0im, e2);               \
      polymod_const_mul(e2, inim[4], p->xpowq4); \
      element_add(e0im, e0im, e2);               \
      element_add(e0im0, e0im0, inim[0]);        \
    } else {                                     \
      polymod_const_mul(e2, inim[1], p->xpowq);  \
      element_neg(e0im, e2);                     \
      polymod_const_mul(e2, inim[2], p->xpowq2); \
      element_sub(e0im, e0im, e2);               \
      polymod_const_mul(e2, inim[3], p->xpowq3); \
      element_sub(e0im, e0im, e2);               \
      polymod_const_mul(e2, inim[4], p->xpowq4); \
      element_sub(e0im, e0im, e2);               \
      element_sub(e0im0, e0im0, inim[0]);        \
    }                                            \
  }
  qpower(1);
  element_set(e3, e0);
  element_set(e0re, element_x(in));
  element_neg(e0im, element_y(in));
  element_mul(e3, e3, e0);
  qpower(-1);
  element_mul(e0, e0, in);
  element_invert(e0, e0);
  element_mul(in, e3, e0);

  element_set(e0, in);
  lucas_even(out, e0, pairing->phikonr);

  element_clear(e0);
  element_clear(e1);
  element_clear(e2);
  element_clear(e3);
  #undef qpower
}

static void (*cc_miller_no_denom_fn)(element_t res, mpz_t q, element_t P,
    element_ptr Qx, element_ptr Qy);

static void cc_pairing(element_ptr out, element_ptr in1, element_ptr in2,
    pairing_t pairing) {
  element_ptr Qbase = in2;
  element_t Qx, Qy;
  mnt_pairing_data_ptr p = pairing->data;

  element_init(Qx, p->Fqd);
  element_init(Qy, p->Fqd);
  //map from twist: (x, y) --> (v^-1 x, v^-(3/2) y)
  //where v is the quadratic nonresidue used to construct the twist
  element_mul(Qx, curve_x_coord(Qbase), p->nqrinv);
  //v^-3/2 = v^-2 * v^1/2
  element_mul(Qy, curve_y_coord(Qbase), p->nqrinv2);
  cc_miller_no_denom_fn(out, pairing->r, in1, Qx, Qy);
  tatepower10(out, out, pairing);
  element_clear(Qx);
  element_clear(Qy);
}

static int cc_is_almost_coddh(element_ptr a, element_ptr b,
    element_ptr c, element_ptr d,
    pairing_t pairing) {
  int res = 0;
  element_t t0, t1, t2;
  element_t cx, cy;
  element_t dx, dy;
  mnt_pairing_data_ptr p = pairing->data;

  element_init(cx, p->Fqd);
  element_init(cy, p->Fqd);
  element_init(dx, p->Fqd);
  element_init(dy, p->Fqd);

  element_init(t0, p->Fqk);
  element_init(t1, p->Fqk);
  element_init(t2, p->Fqk);
  //map from twist: (x, y) --> (v^-1 x, v^-(3/2) y)
  //where v is the quadratic nonresidue used to construct the twist
  element_mul(cx, curve_x_coord(c), p->nqrinv);
  element_mul(dx, curve_x_coord(d), p->nqrinv);
  //v^-3/2 = v^-2 * v^1/2
  element_mul(cy, curve_y_coord(c), p->nqrinv2);
  element_mul(dy, curve_y_coord(d), p->nqrinv2);

  cc_miller_no_denom_fn(t0, pairing->r, a, dx, dy);
  cc_miller_no_denom_fn(t1, pairing->r, b, cx, cy);
  tatepower10(t0, t0, pairing);
  tatepower10(t1, t1, pairing);
  element_mul(t2, t0, t1);
  if (element_is1(t2)) {
    //g, g^x, h, h^-x case
    res = 1;
  } else {
    element_invert(t1, t1);
    element_mul(t2, t0, t1);
    if (element_is1(t2)) {
      //g, g^x, h, h^x case
      res = 1;
    }
  }
  element_clear(cx);
  element_clear(cy);
  element_clear(dx);
  element_clear(dy);
  element_clear(t0);
  element_clear(t1);
  element_clear(t2);
  return res;
}

struct pp_coeff_s {
  element_t a;
  element_t b;
  element_t c;
};
typedef struct pp_coeff_s pp_coeff_t[1];
typedef struct pp_coeff_s *pp_coeff_ptr;

static void g_pairing_pp_init(pairing_pp_t p, element_ptr in1, pairing_t pairing) {
  element_ptr P = in1;
  const element_ptr Px = curve_x_coord(P);
  const element_ptr Py = curve_y_coord(P);
  element_t Z;
  int m;
  mnt_pairing_data_ptr info = pairing->data;
  element_t t0;
  element_t a, b, c;
  field_ptr Fq = info->Fq;
  pp_coeff_t *coeff;
  mpz_ptr q = pairing->r;
  pp_coeff_ptr pp;
  const element_ptr cca = curve_a_coeff(P);
  element_ptr Zx;
  element_ptr Zy;

  #define store_abc() {      \
    element_init(pp->a, Fq); \
    element_init(pp->b, Fq); \
    element_init(pp->c, Fq); \
    element_set(pp->a, a);   \
    element_set(pp->b, b);   \
    element_set(pp->c, c);   \
    pp++;                    \
  }

  //a = -slope_tangent(Z.x, Z.y);
  //b = 1;
  //c = -(Z.y + a * Z.x);
  //but we multiply by 2*Z.y to avoid division

  //a = -Zx * (3 Zx + twicea_2) - a_4;
  //Common curves: a2 = 0 (and cc->a is a_4), so
  //a = -(3 Zx^2 + cc->a)
  //b = 2 * Zy
  //c = -(2 Zy^2 + a Zx);
  #define do_tangent() {    \
    element_square(a, Zx);  \
    element_double(t0, a);  \
    element_add(a, a, t0);  \
    element_add(a, a, cca); \
    element_neg(a, a);      \
                            \
    element_add(b, Zy, Zy); \
                            \
    element_mul(t0, b, Zy); \
    element_mul(c, a, Zx);  \
    element_add(c, c, t0);  \
    element_neg(c, c);      \
                            \
    store_abc();            \
  }

  //a = -(B.y - A.y) / (B.x - A.x);
  //b = 1;
  //c = -(A.y + a * A.x);
  //but we'll multiply by B.x - A.x to avoid division
  #define do_line() {       \
    element_sub(b, Px, Zx); \
    element_sub(a, Zy, Py); \
    element_mul(t0, b, Zy); \
    element_mul(c, a, Zx);  \
    element_add(c, c, t0);  \
    element_neg(c, c);      \
    store_abc();            \
  }

  element_init(Z, P->field);
  element_set(Z, P);
  Zx = curve_x_coord(Z);
  Zy = curve_y_coord(Z);

  element_init(t0, Fq);
  element_init(a, Fq);
  element_init(b, Fq);
  element_init(c, Fq);

  m = mpz_sizeinbase(q, 2) - 2;
  p->data = pbc_malloc(sizeof(pp_coeff_t) * 2 * m);
  coeff = (pp_coeff_t *) p->data;
  pp = coeff[0];

  for(;;) {
    do_tangent();
    if (!m) break;
    element_double(Z, Z);
    if (mpz_tstbit(q, m)) {
      do_line();
      element_add(Z, Z, P);
    }
    m--;
  }

  element_clear(t0);
  element_clear(a);
  element_clear(b);
  element_clear(c);
  element_clear(Z);
  #undef store_abc
  #undef do_tangent
  #undef do_line
}

static void g_pairing_pp_clear(pairing_pp_t p) {
  //TODO: better to store a sentinel value in p->data?
  mpz_ptr q = p->pairing->r;
  int m = mpz_sizeinbase(q, 2) + mpz_popcount(q) - 3;
  int i;
  pp_coeff_t *coeff = (pp_coeff_t *) p->data;
  pp_coeff_ptr pp;
  for (i=0; i<m; i++) {
    pp = coeff[i];
    element_clear(pp->a);
    element_clear(pp->b);
    element_clear(pp->c);
  }
  pbc_free(p->data);
}

static void g_pairing_pp_apply(element_ptr out, element_ptr in2, pairing_pp_t p) {
  mpz_ptr q = p->pairing->r;
  mnt_pairing_data_ptr info = p->pairing->data;
  int m = mpz_sizeinbase(q, 2) - 2;
  pp_coeff_t *coeff = (pp_coeff_t *) p->data;
  pp_coeff_ptr pp = coeff[0];
  element_ptr Qbase = in2;
  element_t e0;
  element_t Qx, Qy;
  element_t v;
  element_init_same_as(e0, out);
  element_init_same_as(v, out);
  element_init(Qx, info->Fqd);
  element_init(Qy, info->Fqd);

  //map from twist: (x, y) --> (v^-1 x, v^-(3/2) y)
  //where v is the quadratic nonresidue used to construct the twist
  element_mul(Qx, curve_x_coord(Qbase), info->nqrinv);
  //v^-3/2 = v^-2 * v^1/2
  element_mul(Qy, curve_y_coord(Qbase), info->nqrinv2);

  element_set1(out);
  for(;;) {
    d_miller_evalfn(e0, pp->a, pp->b, pp->c, Qx, Qy);
    element_mul(out, out, e0);
    pp++;

    if (!m) break;

    if (mpz_tstbit(q, m)) {
      d_miller_evalfn(e0, pp->a, pp->b, pp->c, Qx, Qy);
      element_mul(out, out, e0);
      pp++;
    }
    m--;
    element_square(out, out);
  }
  tatepower10(out, out, p->pairing);

  element_clear(e0);
  element_clear(Qx);
  element_clear(Qy);
  element_clear(v);
}

// in1, in2 are from E(F_q), out from F_q^2
// Compute pairing via elliptic nets (see Stange).
static void g_pairing_ellnet(element_ptr out, element_ptr in1, element_ptr in2,
    pairing_t pairing) {
  mnt_pairing_data_ptr p = pairing->data;

  const element_ptr a = curve_a_coeff(in1);
  const element_ptr b = curve_b_coeff(in1);

  element_ptr x = curve_x_coord(in1);
  element_ptr y = curve_y_coord(in1);

  element_ptr x2 = curve_x_coord(in2);
  element_ptr y2 = curve_y_coord(in2);

  //we map (x2,y2) to (-x2, i y2) before pairing
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
  element_init_same_as(A, out);
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

  //compute A, B, d1

  element_mul(element_x(d0), x2, p->nqrinv);
  element_neg(A, d0);
  element_add(element_item(element_x(A), 0), element_item(element_x(A), 0), x);

  element_double(C, x);
  element_add(element_item(element_x(d0), 0), element_item(element_x(d0), 0), C);

  element_square(dm1, A);
  element_mul(dm1, d0, dm1);

  element_mul(element_y(d1), y2, p->nqrinv2);
  element_set(element_item(element_x(d1), 0), y);

  element_square(d1, d1);
  element_sub(d1, dm1, d1);
  element_invert(B, d1);

  element_invert(A, A);

  element_mul(element_y(d1), y2, p->nqrinv2);
  element_set0(element_x(d1));
  element_neg(element_item(element_x(d1), 0), y);
  element_mul(d1, d1, A);
  element_square(d1, d1);
  element_sub(d1, d0, d1);

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

      polymod_const_mul(element_x(out), t0, element_x(u));
      polymod_const_mul(element_y(out), t0, element_y(u));
      polymod_const_mul(element_x(dm1), s0, element_x(v));
      polymod_const_mul(element_y(dm1), s0, element_y(v));
      element_sub(dm1, dm1, out);

      polymod_const_mul(element_x(out), t1, element_x(u));
      polymod_const_mul(element_y(out), t1, element_y(u));
      polymod_const_mul(element_x(d0), s1, element_x(v));
      polymod_const_mul(element_y(d0), s1, element_y(v));
      element_sub(d0, d0, out);
      element_mul(d0, d0, A);

      polymod_const_mul(element_x(out), t2, element_x(u));
      polymod_const_mul(element_y(out), t2, element_y(u));
      polymod_const_mul(element_x(d1), s2, element_x(v));
      polymod_const_mul(element_y(d1), s2, element_y(v));
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

      polymod_const_mul(element_x(out), tm1, element_x(u));
      polymod_const_mul(element_y(out), tm1, element_y(u));
      polymod_const_mul(element_x(dm1), sm1, element_x(v));
      polymod_const_mul(element_y(dm1), sm1, element_y(v));
      element_sub(dm1, dm1, out);

      polymod_const_mul(element_x(out), t0, element_x(u));
      polymod_const_mul(element_y(out), t0, element_y(u));
      polymod_const_mul(element_x(d0), s0, element_x(v));
      polymod_const_mul(element_y(d0), s0, element_y(v));
      element_sub(d0, d0, out);

      polymod_const_mul(element_x(out), t1, element_x(u));
      polymod_const_mul(element_y(out), t1, element_y(u));
      polymod_const_mul(element_x(d1), s1, element_x(v));
      polymod_const_mul(element_y(d1), s1, element_y(v));
      element_sub(d1, d1, out);
      element_mul(d1, d1, A);
    }
    if (!m) break;
    m--;
  }
  // since c_k lies base field
  // it gets killed by the final powering
  //element_invert(c1, c1);
  //element_mul(element_x(d1), element_x(d1), c1);
  //element_mul(element_y(d1), element_y(d1), c1);

  tatepower10(out, d1, pairing);

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

static void g_pairing_clear(pairing_t pairing) {
  field_clear(pairing->GT);
  mnt_pairing_data_ptr p = pairing->data;

  element_clear(p->xpowq);
  element_clear(p->xpowq2);
  element_clear(p->xpowq3);
  element_clear(p->xpowq4);
  mpz_clear(pairing->phikonr);

  field_clear(p->Etwist);
  field_clear(p->Eq);
  element_clear(p->nqrinv);
  element_clear(p->nqrinv2);
  field_clear(p->Fqk);
  field_clear(p->Fqd);
  field_clear(p->Fqx);
  field_clear(p->Fq);
  field_clear(pairing->Zr);
  mpz_clear(pairing->r);
  pbc_free(p);
}

static void g_pairing_option_set(pairing_t pairing, char *key, char *value) {
  UNUSED_VAR(pairing);
  if (!strcmp(key, "method")) {
    if (!strcmp(value, "miller")) {
      cc_miller_no_denom_fn = cc_miller_no_denom_proj;
    } else if (!strcmp(value, "miller-affine")) {
      cc_miller_no_denom_fn = cc_miller_no_denom_affine;
    } else if (!strcmp(value, "shipsey-stange")) {
      pairing->map = g_pairing_ellnet;
    }
  }
}

static void g_finalpow(element_ptr e) {
  element_t t0;
  element_init_same_as(t0, e->data);
  tatepower10(t0, e->data, e->field->pairing);
  element_set(e->data, t0);
  element_clear(t0);
}

// Computes a curve and sets fp to the field it is defined over using the
// complex multiplication method, where cm holds appropriate data
// (e.g. discriminant, field order).
static void compute_cm_curve(g_param_ptr param, pbc_cm_ptr cm) {
  element_t hp, root;
  field_t fp, fpx;
  field_t cc;

  field_init_fp(fp, cm->q);
  field_init_poly(fpx, fp);
  element_init(hp, fpx);

  mpz_t *coefflist;
  int n = pbc_hilbert(&coefflist, cm->D);

  // Temporarily set the coefficient of x^{n-1} to 1 so hp has degree n - 1,
  // allowing us to use element_item().
  poly_set_coeff1(hp, n - 1);
  int i;
  for (i = 0; i < n; i++) {
    element_set_mpz(element_item(hp, i), coefflist[i]);
  }
  pbc_hilbert_free(coefflist, n);

  //TODO: remove x = 0, 1728 roots
  //TODO: what if there's no roots?
  //printf("hp ");
  //element_out_str(stdout, 0, hp);
  //printf("\n");

  element_init(root, fp);
  poly_findroot(root, hp);
  //printf("root = ");
  //element_out_str(stdout, 0, root);
  //printf("\n");
  element_clear(hp);
  field_clear(fpx);

  //the root is the j-invariant of our desired curve
  field_init_curve_j(cc, root, cm->n, NULL);
  element_clear(root);

  //we may need to twist it however
  {
    // Pick a random point P and twist the curve if it has the wrong order.
    element_t P;
    element_init(P, cc);
    element_random(P);
    element_mul_mpz(P, P, cm->n);
    if (!element_is0(P)) field_reinit_curve_twist(cc);
    element_clear(P);
  }

  mpz_set(param->q, cm->q);
  mpz_set(param->n, cm->n);
  mpz_set(param->h, cm->h);
  mpz_set(param->r, cm->r);
  element_to_mpz(param->a, curve_field_a_coeff(cc));
  element_to_mpz(param->b, curve_field_b_coeff(cc));
  {
    mpz_t z;
    mpz_init(z);
    //compute order of curve in F_q^k
    //n = q - t + 1 hence t = q - n + 1
    mpz_sub(z, param->q, param->n);
    mpz_add_ui(z, z, 1);
    pbc_mpz_trace_n(z, param->q, z, 10);
    mpz_pow_ui(param->nk, param->q, 10);
    mpz_sub_ui(z, z, 1);
    mpz_sub(param->nk, param->nk, z);
    mpz_mul(z, param->r, param->r);
    mpz_divexact(param->hk, param->nk, z);
    mpz_clear(z);
  }
  field_clear(cc);
  field_clear(fp);
}

static void g_init_pairing(pairing_t pairing, void *data) {
  g_param_ptr param = data;
  mnt_pairing_data_ptr p;
  element_t a, b;
  element_t irred;
  int i;

  mpz_init(pairing->r);
  mpz_set(pairing->r, param->r);
  field_init_fp(pairing->Zr, pairing->r);
  pairing->map = cc_pairing;
  pairing->is_almost_coddh = cc_is_almost_coddh;

  p = pairing->data = pbc_malloc(sizeof(mnt_pairing_data_t));
  field_init_fp(p->Fq, param->q);
  element_init(a, p->Fq);
  element_init(b, p->Fq);
  element_set_mpz(a, param->a);
  element_set_mpz(b, param->b);
  field_init_curve_ab(p->Eq, a, b, pairing->r, param->h);

  field_init_poly(p->Fqx, p->Fq);
  element_init(irred, p->Fqx);

  // First set the coefficient of x^5 to 1 so we can call element_item()
  // for the other coefficients.
  poly_set_coeff1(irred, 5);
  for (i=0; i<5; i++) {
    element_set_mpz(element_item(irred, i), param->coeff[i]);
  }

  field_init_polymod(p->Fqd, irred);
  element_clear(irred);

  p->Fqd->nqr = pbc_malloc(sizeof(element_t));
  element_init(p->Fqd->nqr, p->Fqd);
  element_set_mpz(((element_t *) p->Fqd->nqr->data)[0], param->nqr);

  field_init_quadratic(p->Fqk, p->Fqd);

  // Compute phi(k)/r = (q^4 - q^3 + ... + 1)/r.
  {
    element_ptr e = p->xpowq;
    mpz_t z0;
    mpz_ptr q = param->q;
    mpz_ptr z = pairing->phikonr;
    mpz_init(z);
    mpz_init(z0);
    mpz_set_ui(z, 1);
    mpz_sub(z, z, q);
    mpz_mul(z0, q, q);
    mpz_add(z, z, z0);
    mpz_mul(z0, z0, q);
    mpz_sub(z, z, z0);
    mpz_mul(z0, z0, q);
    mpz_add(z, z, z0);
    mpz_clear(z0);
    mpz_divexact(z, z, pairing->r);

    element_init(e, p->Fqd);
    element_init(p->xpowq2, p->Fqd);
    element_init(p->xpowq3, p->Fqd);
    element_init(p->xpowq4, p->Fqd);
    element_set1(((element_t *) e->data)[1]);
    element_pow_mpz(e, e, q);

    element_square(p->xpowq2, p->xpowq);
    element_square(p->xpowq4, p->xpowq2);
    element_mul(p->xpowq3, p->xpowq2, p->xpowq);
  }

  field_init_curve_ab_map(p->Etwist, p->Eq, element_field_to_polymod, p->Fqd, pairing->r, NULL);
  field_reinit_curve_twist(p->Etwist);

  element_init(p->nqrinv, p->Fqd);
  element_invert(p->nqrinv, field_get_nqr(p->Fqd));
  element_init(p->nqrinv2, p->Fqd);
  element_square(p->nqrinv2, p->nqrinv);

  mpz_t ndonr;
  mpz_init(ndonr);
  // ndonr temporarily holds the trace.
  mpz_sub(ndonr, param->q, param->n);
  mpz_add_ui(ndonr, ndonr, 1);
  // Negate because we want the order of the twist.
  mpz_neg(ndonr, ndonr);
  pbc_mpz_curve_order_extn(ndonr, param->q, ndonr, 5);
  mpz_divexact(ndonr, ndonr, param->r);
  field_curve_set_quotient_cmp(p->Etwist, ndonr);
  mpz_clear(ndonr);

  pairing->G1 = p->Eq;
  pairing->G2 = p->Etwist;
  pairing_GT_init(pairing, p->Fqk);
  pairing->finalpow = g_finalpow;

  cc_miller_no_denom_fn = cc_miller_no_denom_affine;
  pairing->option_set = g_pairing_option_set;
  pairing->pp_init = g_pairing_pp_init;
  pairing->pp_clear = g_pairing_pp_clear;
  pairing->pp_apply = g_pairing_pp_apply;

  pairing->clear_func = g_pairing_clear;

  element_clear(a);
  element_clear(b);
}

static void g_init(pbc_param_ptr p) {
  static pbc_param_interface_t interface = {{
    g_clear,
    g_init_pairing,
    g_out_str,
  }};
  p->api = interface;
  g_param_ptr param = p->data = pbc_malloc(sizeof(*param));
  mpz_init(param->q);
  mpz_init(param->n);
  mpz_init(param->h);
  mpz_init(param->r);
  mpz_init(param->a);
  mpz_init(param->b);
  mpz_init(param->nk);
  mpz_init(param->hk);
  param->coeff = NULL;
  mpz_init(param->nqr);
}

// Public interface:

int pbc_param_init_g(pbc_param_ptr par, struct symtab_s *tab) {
  g_init(par);
  g_param_ptr p = par->data;
  char s[80];

  int err = 0;
  err += lookup_mpz(p->q, tab, "q");
  err += lookup_mpz(p->n, tab, "n");
  err += lookup_mpz(p->h, tab, "h");
  err += lookup_mpz(p->r, tab, "r");
  err += lookup_mpz(p->a, tab, "a");
  err += lookup_mpz(p->b, tab, "b");
  err += lookup_mpz(p->nk, tab, "nk");
  err += lookup_mpz(p->hk, tab, "hk");
  err += lookup_mpz(p->nqr, tab, "nqr");

  p->coeff = pbc_realloc(p->coeff, sizeof(mpz_t) * 5);
  int i;
  for (i = 0; i < 5; i++) {
    sprintf(s, "coeff%d", i);
    mpz_init(p->coeff[i]);
    err += lookup_mpz(p->coeff[i], tab, s);
  }
  return err;
}

void pbc_param_init_g_gen(pbc_param_t p, pbc_cm_ptr cm) {
  g_init(p);
  g_param_ptr param = p->data;
  field_t Fq, Fqx, Fqd;
  element_t irred, nqr;
  int i;

  compute_cm_curve(param, cm);

  field_init_fp(Fq, param->q);
  field_init_poly(Fqx, Fq);
  element_init(irred, Fqx);
  do {
    poly_random_monic(irred, 5);
  } while (!poly_is_irred(irred));
  field_init_polymod(Fqd, irred);

  // Find a quadratic nonresidue of Fqd lying in Fq.
  element_init(nqr, Fqd);
  do {
    element_random(((element_t *) nqr->data)[0]);
  } while (element_is_sqr(nqr));

  param->coeff = pbc_realloc(param->coeff, sizeof(mpz_t) * 5);

  for (i=0; i<5; i++) {
    mpz_init(param->coeff[i]);
    element_to_mpz(param->coeff[i], element_item(irred, i));
  }
  element_to_mpz(param->nqr, ((element_t *) nqr->data)[0]);

  element_clear(nqr);
  element_clear(irred);

  field_clear(Fqx);
  field_clear(Fqd);
  field_clear(Fq);
}
