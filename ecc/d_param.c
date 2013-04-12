// Type D pairings, aka MNT curves.

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
#include "pbc_d_param.h"
#include "ecc/param.h"

struct d_param_s {
  mpz_t q;       // curve defined over F_q
  mpz_t n;       // has order n (= q - t + 1) in F_q
  mpz_t h;       // h * r = n, r is prime
  mpz_t r;
  mpz_t a, b;    // curve equation is y^2 = x^3 + ax + b
  int k;         // embedding degree
  mpz_t nk;      // order of curve over F_q^k
  mpz_t hk;      // hk * r^2 = nk
  mpz_t *coeff;  // coefficients of polynomial used to extend F_q by k/2
  mpz_t nqr;     // a quadratic nonresidue in F_q^d that lies in F_q
};

typedef struct d_param_s d_param_t[1];
typedef struct d_param_s *d_param_ptr;

// Per-pairing data.
typedef struct {
  field_t Fq, Fqx, Fqd, Fqk;  // The fields F_q, F_q[x], F_q^d, F_q^k.
  field_t Eq, Etwist;         // The curves E(F_q) and E'(F_q^d).
  // Let v be the quadratic nonresidue used to construct F_q^k from F_q^d,
  // namely Fqk = Fqd[sqrt(v)].
  element_t nqrinv, nqrinv2;  // The constants v^-1 and v^-2.
  mpz_t tateexp;              // The Tate exponent,
                              // to standardize coset representatives.
  int k;                      // The embedding degree, usually 6.
  // Let x be the element used to build Fqd from Fq, i.e. Fqd = Fq[x].
  element_t xpowq, xpowq2;    // x^q and x^{2q} in F_q^d.
} *pptr;

static void d_clear(void *data) {
  d_param_ptr param = data;
  int d = param->k / 2;
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
  for (i=0; i<d; i++) {
    mpz_clear(param->coeff[i]);
  }
  pbc_free(param->coeff);
  pbc_free(data);
}

static void d_out_str(FILE *stream, void *data) {
  d_param_ptr p = data;
  int d = p->k / 2;
  int i;
  char s[8];
  param_out_type(stream, "d");
  param_out_mpz(stream, "q", p->q);
  param_out_mpz(stream, "n", p->n);
  param_out_mpz(stream, "h", p->h);
  param_out_mpz(stream, "r", p->r);
  param_out_mpz(stream, "a", p->a);
  param_out_mpz(stream, "b", p->b);
  param_out_int(stream, "k", p->k);
  param_out_mpz(stream, "nk", p->nk);
  param_out_mpz(stream, "hk", p->hk);
  for (i=0; i<d; i++) {
    sprintf(s, "coeff%d", i);
    param_out_mpz(stream, s, p->coeff[i]);
  }
  param_out_mpz(stream, "nqr", p->nqr);
}

// Define l = aX + bY + c where a, b, c are in Fq.
// Compute e0 = l(Q) specialized for the case when Q has the form
// (Qx, Qy * sqrt(v)) where Qx, Qy are in Fqd and v is the quadratic nonresidue
// used to construct the quadratic field extension Fqk of Fqd.
static inline void d_miller_evalfn(element_t e0,
    element_t a, element_t b, element_t c, element_t Qx, element_t Qy) {
  element_ptr re_out = element_x(e0);
  element_ptr im_out = element_y(e0);

  int i;
  int d = polymod_field_degree(re_out->field);
  for (i = 0; i < d; i++) {
    element_mul(element_item(re_out, i), element_item(Qx, i), a);
    element_mul(element_item(im_out, i), element_item(Qy, i), b);
  }
  element_add(element_item(re_out, 0), element_item(re_out, 0), c);
}

// Miller's algorithm, assuming we can ignore the denominator. We can do this
// with careful group selection when the embedding degree is even. See thesis.
// This version uses projective coordinates, which don't seem much faster.
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
    /* t4 holds e7 t3^2. */                     \
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

// Same as above, but with affine coordinates.
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
  void do_vertical() {
    mapbase(e0, Z->x);
    element_sub(e0, Qx, e0);
    element_mul(v, v, e0);
  }
  */

  #define do_tangent() {                  \
    /* a = -(3 Zx^2 + cc->a) */           \
    /* b = 2 * Zy */                      \
    /* c = -(2 Zy^2 + a Zx); */           \
                                          \
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

  #define do_line() {                                     \
    /* a = -(B.y - A.y) / (B.x - A.x); */                 \
    /* b = 1; */                                          \
    /* c = -(A.y + a * A.x); */                           \
    /* but we multiply by B.x - A.x to avoid division. */ \
                                                          \
    element_sub(b, Px, Zx);                               \
    element_sub(a, Zy, Py);                               \
    element_mul(t0, b, Zy);                               \
    element_mul(c, a, Zx);                                \
    element_add(c, c, t0);                                \
    element_neg(c, c);                                    \
                                                          \
    d_miller_evalfn(e0, a, b, c, Qx, Qy);                 \
    element_mul(v, v, e0);                                \
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

static void (*cc_miller_no_denom_fn)(element_t res, mpz_t q, element_t P,
    element_ptr Qx, element_ptr Qy);

static void d_pairing_option_set(pairing_t pairing, char *key, char *value) {
  UNUSED_VAR(pairing);
  if (!strcmp(key, "method")) {
    if (!strcmp(value, "miller")) {
      cc_miller_no_denom_fn = cc_miller_no_denom_proj;
    } else if (!strcmp(value, "miller-affine")) {
      cc_miller_no_denom_fn = cc_miller_no_denom_affine;
    }
  }
}

// Requires cofactor is even. TODO: This seems to contradict a comment below.
// Requires in != out.
// Mangles in.
static void lucas_even(element_ptr out, element_ptr in, mpz_t cofactor) {
  if (element_is1(in)) {
    element_set(out, in);
    return;
  }
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

  // Assume cofactor = (q^2 - q + 1) / r is odd
  // thus v1 = V_k, v0 = V_{k-1}
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

// The final powering, where we standardize the coset representative.
static void cc_tatepower(element_ptr out, element_ptr in, pairing_t pairing) {
  pptr p = pairing->data;
  #define qpower(sign) {                         \
    polymod_const_mul(e2, inre[1], p->xpowq);    \
    element_set(e0re, e2);                       \
    polymod_const_mul(e2, inre[2], p->xpowq2);   \
    element_add(e0re, e0re, e2);                 \
    element_add(e0re0, e0re0, inre[0]);          \
                                                 \
    if (sign > 0) {                              \
      polymod_const_mul(e2, inim[1], p->xpowq);  \
      element_set(e0im, e2);                     \
      polymod_const_mul(e2, inim[2], p->xpowq2); \
      element_add(e0im, e0im, e2);               \
      element_add(e0im0, e0im0, inim[0]);        \
    } else {                                     \
      polymod_const_mul(e2, inim[1], p->xpowq);  \
      element_neg(e0im, e2);                     \
      polymod_const_mul(e2, inim[2], p->xpowq2); \
      element_sub(e0im, e0im, e2);               \
      element_sub(e0im0, e0im0, inim[0]);        \
    }                                            \
  }
  if (p->k == 6) {
    // See thesis, section 6.9, "The Final Powering", which gives a formula
    // for the first step of the final powering when Fq6 has been implemented
    // as a quadratic extension on top of a cubic extension.
    element_t e0, e2, e3;
    element_init(e0, p->Fqk);
    element_init(e2, p->Fqd);
    element_init(e3, p->Fqk);
    element_ptr e0re = element_x(e0);
    element_ptr e0im = element_y(e0);
    element_ptr e0re0 = ((element_t *) e0re->data)[0];
    element_ptr e0im0 = ((element_t *) e0im->data)[0];
    element_t *inre = element_x(in)->data;
    element_t *inim = element_y(in)->data;
    // Expressions in the formula are similar, hence the following function.
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
    // We use Lucas sequences to complete the final powering.
    lucas_even(out, e0, pairing->phikonr);

    element_clear(e0);
    element_clear(e2);
    element_clear(e3);
  } else {
    element_pow_mpz(out, in, p->tateexp);
  }
  #undef qpower
}

static void cc_finalpow(element_t e) {
  cc_tatepower(e->data, e->data, e->field->pairing);
}

static void cc_pairing(element_ptr out, element_ptr in1, element_ptr in2,
    pairing_t pairing) {
  element_ptr Qbase = in2;
  element_t Qx, Qy;
  pptr p = pairing->data;

  element_init(Qx, p->Fqd);
  element_init(Qy, p->Fqd);
  // Twist: (x, y) --> (v^-1 x, v^-(3/2) y)
  // where v is the quadratic nonresidue used to construct the twist.
  element_mul(Qx, curve_x_coord(Qbase), p->nqrinv);
  // v^-3/2 = v^-2 * v^1/2
  element_mul(Qy, curve_y_coord(Qbase), p->nqrinv2);
  cc_miller_no_denom_fn(out, pairing->r, in1, Qx, Qy);
  cc_tatepower(out, out, pairing);
  element_clear(Qx);
  element_clear(Qy);
}


//do many millers at one time with affine coordinates.
static void cc_millers_no_denom_affine(element_t res, mpz_t q, element_t P[],
    element_t Qx[], element_t Qy[], int n_prod) {
  int m, i;
  element_t v;
  element_t a, b, c;
  element_t t0;
  element_t e0;
  const element_ptr cca = curve_a_coeff(P[0]);
  element_ptr Px, Py;
  element_t* Z = pbc_malloc(sizeof(element_t)*n_prod);
  element_ptr Zx, Zy;

  /* TODO: when exactly is this not needed?
  void do_vertical() {
    mapbase(e0, Z->x);
    element_sub(e0, Qx, e0);
    element_mul(v, v, e0);
  }
  */

  #define do_tangents() {                         \
    /* a = -(3 Zx^2 + cc->a) */                   \
    /* b = 2 * Zy */                              \
    /* c = -(2 Zy^2 + a Zx); */                   \
    for(i=0; i<n_prod; i++){                      \
      Px = curve_x_coord(P[i]);                   \
      Py = curve_y_coord(P[i]);                   \
      Zx = curve_x_coord(Z[i]);                   \
      Zy = curve_y_coord(Z[i]);                   \
                                                  \
      element_square(a, Zx);                      \
      element_mul_si(a, a, 3);                    \
      element_add(a, a, cca);                     \
      element_neg(a, a);                          \
                                                  \
      element_add(b, Zy, Zy);                     \
                                                  \
      element_mul(t0, b, Zy);                     \
      element_mul(c, a, Zx);                      \
      element_add(c, c, t0);                      \
      element_neg(c, c);                          \
                                                  \
      d_miller_evalfn(e0, a, b, c, Qx[i], Qy[i]); \
      element_mul(v, v, e0);                      \
    }                                             \
  }

  #define do_lines() {                                    \
    /* a = -(B.y - A.y) / (B.x - A.x); */                 \
    /* b = 1; */                                          \
    /* c = -(A.y + a * A.x); */                           \
    /* but we multiply by B.x - A.x to avoid division. */ \
    for(i=0; i<n_prod; i++){                              \
      Px = curve_x_coord(P[i]);                           \
      Py = curve_y_coord(P[i]);                           \
      Zx = curve_x_coord(Z[i]);                           \
      Zy = curve_y_coord(Z[i]);                           \
                                                          \
      element_sub(b, Px, Zx);                             \
      element_sub(a, Zy, Py);                             \
      element_mul(t0, b, Zy);                             \
      element_mul(c, a, Zx);                              \
      element_add(c, c, t0);                              \
      element_neg(c, c);                                  \
                                                          \
      d_miller_evalfn(e0, a, b, c, Qx[i], Qy[i]);         \
      element_mul(v, v, e0);                              \
    }                                                     \
  }

  Px= curve_x_coord(P[0]); //temporally used to initial a,b, c and etc.
  element_init(a, Px->field);
  element_init(b, a->field);
  element_init(c, a->field);
  element_init(t0, a->field);
  element_init(e0, res->field);

  element_init(v, res->field);
  for(i=0; i<n_prod; i++){
    element_init(Z[i], P[i]->field);
    element_set(Z[i], P[i]);
  }

  element_set1(v);
  m = mpz_sizeinbase(q, 2) - 2;

  for(;;) {
    do_tangents();

    if (!m) break;
    element_multi_double(Z, Z, n_prod); //Z_i=Z_i+Z_i for all i.

    if (mpz_tstbit(q, m)) {
      do_lines();
      element_multi_add(Z, Z, P, n_prod); //Z_i=Z_i+P_i for all i.
    }
    m--;
    element_square(v, v);
  }

  element_set(res, v);

  element_clear(v);
  for(i=0; i<n_prod; i++){
    element_clear(Z[i]);
  }
  pbc_free(Z);
  element_clear(a);
  element_clear(b);
  element_clear(c);
  element_clear(t0);
  element_clear(e0);
  #undef do_tangents
  #undef do_lines
}


void cc_pairings_affine(element_ptr out, element_t in1[], element_t in2[],
        int n_prod, pairing_t pairing) {
  element_ptr Qbase;
  element_t* Qx = pbc_malloc(sizeof(element_t)*n_prod);
  element_t* Qy = pbc_malloc(sizeof(element_t)*n_prod);
  pptr p = pairing->data;
  int i;
  for(i=0; i<n_prod; i++){
          element_init(Qx[i], p->Fqd);
          element_init(Qy[i], p->Fqd);
        Qbase = in2[i];
          // Twist: (x, y) --> (v^-1 x, v^-(3/2) y)
          // where v is the quadratic nonresidue used to construct the twist.
          element_mul(Qx[i], curve_x_coord(Qbase), p->nqrinv);
          // v^-3/2 = v^-2 * v^1/2
          element_mul(Qy[i], curve_y_coord(Qbase), p->nqrinv2);
  }
  cc_millers_no_denom_affine(out, pairing->r, in1, Qx, Qy, n_prod);
  cc_tatepower(out, out, pairing);

  for(i=0; i<n_prod; i++){
          element_clear(Qx[i]);
                element_clear(Qy[i]);
  }
  pbc_free(Qx);
  pbc_free(Qy);
}


static int cc_is_almost_coddh(element_ptr a, element_ptr b,
    element_ptr c, element_ptr d,
    pairing_t pairing) {
  int res = 0;
  element_t t0, t1, t2;
  element_t cx, cy;
  element_t dx, dy;
  pptr p = pairing->data;

  element_init(cx, p->Fqd);
  element_init(cy, p->Fqd);
  element_init(dx, p->Fqd);
  element_init(dy, p->Fqd);

  element_init(t0, p->Fqk);
  element_init(t1, p->Fqk);
  element_init(t2, p->Fqk);
  // Twist: (x, y) --> (v^-1 x, v^-(3/2) y)
  // where v is the quadratic nonresidue used to construct the twist.
  element_mul(cx, curve_x_coord(c), p->nqrinv);
  element_mul(dx, curve_x_coord(d), p->nqrinv);
  // v^-3/2 = v^-2 * v^1/2
  element_mul(cy, curve_y_coord(c), p->nqrinv2);
  element_mul(dy, curve_y_coord(d), p->nqrinv2);

  cc_miller_no_denom_fn(t0, pairing->r, a, dx, dy);
  cc_miller_no_denom_fn(t1, pairing->r, b, cx, cy);
  cc_tatepower(t0, t0, pairing);
  cc_tatepower(t1, t1, pairing);
  element_mul(t2, t0, t1);
  if (element_is1(t2)) res = 1;    // We were given g, g^x, h, h^-x.
  else {
    // Cheaply check the other case.
    element_invert(t1, t1);
    element_mul(t2, t0, t1);
    if (element_is1(t2)) res = 1;  // We were given g, g^x, h, h^x.
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

static void d_pairing_pp_init(pairing_pp_t p, element_ptr in1, pairing_t pairing) {
  element_ptr P = in1;
  const element_ptr Px = curve_x_coord(P);
  const element_ptr Py = curve_y_coord(P);
  element_t Z;
  int m;
  pptr info = pairing->data;
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

  #define do_tangent() {                               \
    /* a = -slope_tangent(Z.x, Z.y); */                \
    /* b = 1; */                                       \
    /* c = -(Z.y + a * Z.x); */                        \
    /* but we multiply by 2*Z.y to avoid division. */  \
                                                       \
    /* a = -Zx * (3 Zx + twicea_2) - a_4; */           \
    /* Common curves: a2 = 0 (and cc->a is a_4), so */ \
    /* a = -(3 Zx^2 + cc->a) */                        \
    /* b = 2 * Zy */                                   \
    /* c = -(2 Zy^2 + a Zx); */                        \
                                                       \
    element_square(a, Zx);                             \
    element_double(t0, a);                             \
    element_add(a, a, t0);                             \
    element_add(a, a, cca);                            \
    element_neg(a, a);                                 \
                                                       \
    element_add(b, Zy, Zy);                            \
                                                       \
    element_mul(t0, b, Zy);                            \
    element_mul(c, a, Zx);                             \
    element_add(c, c, t0);                             \
    element_neg(c, c);                                 \
                                                       \
    store_abc();                                       \
  }

  #define do_line() {                                       \
    /* a = -(B.y - A.y) / (B.x - A.x); */                   \
    /* b = 1; */                                            \
    /* c = -(A.y + a * A.x); */                             \
    /* but we'll multiply by B.x - A.x to avoid division */ \
                                                            \
    element_sub(b, Px, Zx);                                 \
    element_sub(a, Zy, Py);                                 \
    element_mul(t0, b, Zy);                                 \
    element_mul(c, a, Zx);                                  \
    element_add(c, c, t0);                                  \
    element_neg(c, c);                                      \
                                                            \
    store_abc();                                            \
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

static void d_pairing_pp_clear(pairing_pp_t p) {
  // TODO: Better to store a sentinel value in p->data?
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

static void d_pairing_pp_apply(element_ptr out, element_ptr in2,
    pairing_pp_t p) {
  mpz_ptr q = p->pairing->r;
  pptr info = p->pairing->data;
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

  // Twist: (x, y) --> (v^-1 x, v^-(3/2) y)
  // where v is the quadratic nonresidue used to construct the twist
  element_mul(Qx, curve_x_coord(Qbase), info->nqrinv);
  // v^-3/2 = v^-2 * v^1/2
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
  cc_tatepower(out, out, p->pairing);

  element_clear(e0);
  element_clear(Qx);
  element_clear(Qy);
  element_clear(v);
}

static void d_pairing_clear(pairing_t pairing) {
  field_clear(pairing->GT);
  pptr p = pairing->data;

  if (p->k == 6) {
    element_clear(p->xpowq);
    element_clear(p->xpowq2);
    mpz_clear(pairing->phikonr);
  } else {
    mpz_clear(p->tateexp);
  }

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

static void d_init_pairing(pairing_ptr pairing, void *data) {
  d_param_ptr param = data;
  pptr p;
  element_t a, b;
  element_t irred;
  int d = param->k / 2;
  int i;

  if (param->k % 2) pbc_die("k must be even");

  mpz_init(pairing->r);
  mpz_set(pairing->r, param->r);
  field_init_fp(pairing->Zr, pairing->r);
  pairing->map = cc_pairing;
  pairing->prod_pairings = cc_pairings_affine;
  pairing->is_almost_coddh = cc_is_almost_coddh;

  p = pairing->data = pbc_malloc(sizeof(*p));
  field_init_fp(p->Fq, param->q);
  element_init(a, p->Fq);
  element_init(b, p->Fq);
  element_set_mpz(a, param->a);
  element_set_mpz(b, param->b);
  field_init_curve_ab(p->Eq, a, b, pairing->r, param->h);

  field_init_poly(p->Fqx, p->Fq);
  element_init(irred, p->Fqx);
  poly_set_coeff1(irred, d);
  for (i = 0; i < d; i++) {
    element_set_mpz(element_item(irred, i), param->coeff[i]);
  }

  field_init_polymod(p->Fqd, irred);
  element_clear(irred);

  p->Fqd->nqr = pbc_malloc(sizeof(element_t));
  element_init(p->Fqd->nqr, p->Fqd);
  element_set_mpz(((element_t *) p->Fqd->nqr->data)[0], param->nqr);

  field_init_quadratic(p->Fqk, p->Fqd);

  // Compute constants involved in the final powering.
  if (param->k == 6) {
    mpz_ptr q = param->q;
    mpz_ptr z = pairing->phikonr;
    mpz_init(z);
    mpz_mul(z, q, q);
    mpz_sub(z, z, q);
    mpz_add_ui(z, z, 1);
    mpz_divexact(z, z, pairing->r);

    element_ptr e = p->xpowq;
    element_init(e, p->Fqd);
    element_set1(((element_t *) e->data)[1]);
    element_pow_mpz(e, e, q);

    element_init(p->xpowq2, p->Fqd);
    element_square(p->xpowq2, e);
  } else {
    mpz_init(p->tateexp);
    mpz_sub_ui(p->tateexp, p->Fqk->order, 1);
    mpz_divexact(p->tateexp, p->tateexp, pairing->r);
  }

  field_init_curve_ab_map(p->Etwist, p->Eq, element_field_to_polymod, p->Fqd, pairing->r, NULL);
  field_reinit_curve_twist(p->Etwist);

  mpz_t ndonr;
  mpz_init(ndonr);
  // ndonr temporarily holds the trace.
  mpz_sub(ndonr, param->q, param->n);
  mpz_add_ui(ndonr, ndonr, 1);
  // Negate it because we want the trace of the twist.
  mpz_neg(ndonr, ndonr);
  pbc_mpz_curve_order_extn(ndonr, param->q, ndonr, d);
  mpz_divexact(ndonr, ndonr, param->r);
  field_curve_set_quotient_cmp(p->Etwist, ndonr);
  mpz_clear(ndonr);

  element_init(p->nqrinv, p->Fqd);
  element_invert(p->nqrinv, field_get_nqr(p->Fqd));
  element_init(p->nqrinv2, p->Fqd);
  element_square(p->nqrinv2, p->nqrinv);

  pairing->G1 = p->Eq;
  pairing->G2 = p->Etwist;

  p->k = param->k;
  pairing_GT_init(pairing, p->Fqk);
  pairing->finalpow = cc_finalpow;

  // By default use affine coordinates.
  cc_miller_no_denom_fn = cc_miller_no_denom_affine;
  pairing->option_set = d_pairing_option_set;
  pairing->pp_init = d_pairing_pp_init;
  pairing->pp_clear = d_pairing_pp_clear;
  pairing->pp_apply = d_pairing_pp_apply;

  pairing->clear_func = d_pairing_clear;

  element_clear(a);
  element_clear(b);
}

// Computes a curve and sets fp to the field it is defined over using the
// complex multiplication method, where cm holds the appropriate information
// (e.g. discriminant, field order).
static void compute_cm_curve(d_param_ptr param, pbc_cm_ptr cm) {
  element_t hp, root;
  field_t fp, fpx;
  field_t cc;

  field_init_fp(fp, cm->q);
  field_init_poly(fpx, fp);
  element_init(hp, fpx);

  mpz_t *coefflist;
  int n = pbc_hilbert(&coefflist, cm->D);

  // Temporarily set the coefficient of x^{n-1} to 1 so hp has degree n - 1,
  // allowing us to use poly_coeff().
  poly_set_coeff1(hp, n - 1);
  int i;
  for (i = 0; i < n; i++) {
    element_set_mpz(element_item(hp, i), coefflist[i]);
  }
  pbc_hilbert_free(coefflist, n);

  // TODO: Remove x = 0, 1728 roots.
  // TODO: What if there are no roots?
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

  // The root is the j-invariant of the desired curve.
  field_init_curve_j(cc, root, cm->n, NULL);
  element_clear(root);

  // We may need to twist it.
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
  param->k = cm->k;
  {
    mpz_t z;
    mpz_init(z);
    // Compute order of curve in F_q^k.
    // n = q - t + 1 hence t = q - n + 1
    mpz_sub(z, param->q, param->n);
    mpz_add_ui(z, z, 1);
    pbc_mpz_trace_n(z, param->q, z, param->k);
    mpz_pow_ui(param->nk, param->q, param->k);
    mpz_sub_ui(z, z, 1);
    mpz_sub(param->nk, param->nk, z);
    mpz_mul(z, param->r, param->r);
    mpz_divexact(param->hk, param->nk, z);
    mpz_clear(z);
  }
  field_clear(cc);
  field_clear(fp);
}

static void d_param_init(pbc_param_ptr p) {
  static pbc_param_interface_t interface = {{
    d_clear,
    d_init_pairing,
    d_out_str,
  }};
  p->api = interface;
  d_param_ptr param = p->data = pbc_malloc(sizeof(*param));
  mpz_init(param->q);
  mpz_init(param->n);
  mpz_init(param->h);
  mpz_init(param->r);
  mpz_init(param->a);
  mpz_init(param->b);
  mpz_init(param->nk);
  mpz_init(param->hk);
  param->k = 0;
  param->coeff = NULL;
  mpz_init(param->nqr);
}

// Public interface:

int pbc_param_init_d(pbc_param_ptr par, struct symtab_s *tab) {
  d_param_init(par);
  d_param_ptr p = par->data;
  char s[80];
  int i, d;

  int err = 0;
  err += lookup_mpz(p->q, tab, "q");
  err += lookup_mpz(p->n, tab, "n");
  err += lookup_mpz(p->h, tab, "h");
  err += lookup_mpz(p->r, tab, "r");
  err += lookup_mpz(p->a, tab, "a");
  err += lookup_mpz(p->b, tab, "b");
  err += lookup_int(&p->k, tab, "k");
  err += lookup_mpz(p->nk, tab, "nk");
  err += lookup_mpz(p->hk, tab, "hk");
  err += lookup_mpz(p->nqr, tab, "nqr");

  d = p->k / 2;
  p->coeff = pbc_realloc(p->coeff, sizeof(mpz_t) * d);
  for (i=0; i<d; i++) {
    sprintf(s, "coeff%d", i);
    mpz_init(p->coeff[i]);
    err += lookup_mpz(p->coeff[i], tab, s);
  }
  return err;
}

void pbc_param_init_d_gen(pbc_param_ptr p, pbc_cm_ptr cm) {
  d_param_init(p);
  d_param_ptr param = p->data;
  field_t Fq, Fqx, Fqd;
  element_t irred, nqr;
  int d = cm->k / 2;
  int i;

  compute_cm_curve(param, cm);

  field_init_fp(Fq, param->q);
  field_init_poly(Fqx, Fq);
  element_init(irred, Fqx);
  do {
    poly_random_monic(irred, d);
  } while (!poly_is_irred(irred));
  field_init_polymod(Fqd, irred);

  // Find a quadratic nonresidue of Fqd lying in Fq.
  element_init(nqr, Fqd);
  do {
    element_random(((element_t *) nqr->data)[0]);
  } while (element_is_sqr(nqr));

  param->coeff = pbc_realloc(param->coeff, sizeof(mpz_t) * d);

  for (i=0; i<d; i++) {
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
