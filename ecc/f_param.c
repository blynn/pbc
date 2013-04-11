#include <stdarg.h>
#include <stdio.h>
#include <stdint.h> // for intptr_t
#include <stdlib.h>
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_field.h"
#include "pbc_fp.h"
#include "pbc_fieldquadratic.h"
#include "pbc_param.h"
#include "pbc_pairing.h"
#include "pbc_poly.h"
#include "pbc_curve.h"
#include "pbc_memory.h"
#include "pbc_f_param.h"
#include "ecc/param.h"

struct f_param_s {
    mpz_t q; // Curve defined over F_q.
    mpz_t r; // The order of the curve.
    mpz_t b; // E: y^2 = x^3 + b
    mpz_t beta; //beta is a quadratic nonresidue in Fq
        //we use F_q^2 = F_q[sqrt(beta)]
    mpz_t alpha0, alpha1;
        //the polynomial x^6 + alpha0 + alpha1 sqrt(beta)
        //is irreducible over F_q^2[x], so
        //we can extend F_q^2 to F_q^12 using the
        //sixth root of -(alpha0 + alpha1 sqrt(beta))
};
typedef struct f_param_s f_param_t[1];
typedef struct f_param_s *f_param_ptr;

// TODO: we never use phikonr so don't bother computing it,
// but one day other routines might need it
struct f_pairing_data_s {
  field_t Fq, Fq2, Fq2x, Fq12;
  field_t Eq, Etwist;
  element_t negalpha;
  element_t negalphainv;
  mpz_t tateexp;

  //for tate exponentiation speedup:
  //x^{q^k} for various k
  element_t xpowq2, xpowq6, xpowq8;
};
typedef struct f_pairing_data_s f_pairing_data_t[1];
typedef struct f_pairing_data_s *f_pairing_data_ptr;

static void f_clear(void *data) {
  f_param_ptr fp = data;
  mpz_clear(fp->q);
  mpz_clear(fp->r);
  mpz_clear(fp->b);
  mpz_clear(fp->beta);
  mpz_clear(fp->alpha0);
  mpz_clear(fp->alpha1);
  pbc_free(data);
}

static void f_out_str(FILE *stream, void *data) {
  f_param_ptr p = data;
  param_out_type(stream, "f");
  param_out_mpz(stream, "q", p->q);
  param_out_mpz(stream, "r", p->r);
  param_out_mpz(stream, "b", p->b);
  param_out_mpz(stream, "beta", p->beta);
  param_out_mpz(stream, "alpha0", p->alpha0);
  param_out_mpz(stream, "alpha1", p->alpha1);
}

static void tryminusx(mpz_ptr q, mpz_ptr x) {
  //36x4 - 36x3 + 24x2 - 6x + 1
  //= ((36(x - 1)x + 24)x - 6)x + 1
  mpz_sub_ui(q, x, 1);
  mpz_mul(q, q, x);
  mpz_mul_ui(q, q, 36);
  mpz_add_ui(q, q, 24);
  mpz_mul(q, q, x);
  mpz_sub_ui(q, q, 6);
  mpz_mul(q, q, x);
  mpz_add_ui(q, q, 1);
}

static void tryplusx(mpz_ptr q, mpz_ptr x) {
  //36x4 + 36x3 + 24x2 + 6x + 1
  //= ((36(x + 1)x + 24)x + 6)x + 1
  mpz_add_ui(q, x, 1);
  mpz_mul(q, q, x);
  mpz_mul_ui(q, q, 36);
  mpz_add_ui(q, q, 24);
  mpz_mul(q, q, x);
  mpz_add_ui(q, q, 6);
  mpz_mul(q, q, x);
  mpz_add_ui(q, q, 1);
}

static void cc_miller_no_denom(element_t res, mpz_t q, element_t P,
    element_ptr Qx, element_ptr Qy, element_t negalpha) {
  int m;
  element_t v;
  element_t Z;
  element_t a, b, c;
  element_t t0;
  element_t e0, e1;
  element_ptr Zx, Zy;
  const element_ptr Px = curve_x_coord(P);
  const element_ptr Py = curve_y_coord(P);

  #define do_term(i, j, k, flag) {                                \
    element_ptr e2;                                               \
    e2 = element_item(e0, i);                                     \
    element_mul(e1, element_item(v, j), Qx);                      \
    if (flag == 1) element_mul(e1, e1, negalpha);                 \
    element_mul(element_x(e1), element_x(e1), a);                 \
    element_mul(element_y(e1), element_y(e1), a);                 \
    element_mul(e2, element_item(v, k), Qy);                      \
    element_mul(element_x(e2), element_x(e2), b);                 \
    element_mul(element_y(e2), element_y(e2), b);                 \
    element_add(e2, e2, e1);                                      \
    if (flag == 2) element_mul(e2, e2, negalpha);                 \
    element_mul(element_x(e1), element_x(element_item(v, i)), c); \
    element_mul(element_y(e1), element_y(element_item(v, i)), c); \
    element_add(e2, e2, e1);                                      \
  }

  // a, b, c lie in Fq
  // Qx, Qy lie in Fq^2
  // Qx is coefficient of x^4
  // Qy is coefficient of x^3
  //
  // computes v *= (a Qx x^4 + b Qy x^3 + c)
  //
  // recall x^6 = -alpha thus
  // x^4 (u0 + u1 x^1 + ... + u5 x^5) =
  // u0 x^4 + u1 x^5
  // - alpha u2 - alpha u3 x - alpha u4 x^2 - alpha u5 x^3
  // and
  // x^4 (u0 + u1 x^1 + ... + u5 x^5) =
  // u0 x^3 + u1 x^4 + u2 x^5
  // - alpha u3 - alpha u4 x - alpha u5 x^2
  #define f_miller_evalfn() { \
    do_term(0, 2, 3, 2);      \
    do_term(1, 3, 4, 2);      \
    do_term(2, 4, 5, 2);      \
    do_term(3, 5, 0, 1);      \
    do_term(4, 0, 1, 0);      \
    do_term(5, 1, 2, 0);      \
    element_set(v, e0);       \
  }
  /*
    element_ptr e1;

    e1 = element_item(e0, 4);

    element_mul(element_x(e1), element_x(Qx), a);
    element_mul(element_y(e1), element_y(Qx), a);

    e1 = element_item(e0, 3);

    element_mul(element_x(e1), element_x(Qy), b);
    element_mul(element_y(e1), element_y(Qy), b);

    element_set(element_x(element_item(e0, 0)), c);

    element_mul(v, v, e0);
   */

  //a = -3 Zx^2 since cc->a is 0 for D = 3
  //b = 2 * Zy
  //c = -(2 Zy^2 + a Zx);
  #define do_tangent() {     \
    element_square(a, Zx);   \
    element_mul_si(a, a, 3); \
    element_neg(a, a);       \
                             \
    element_add(b, Zy, Zy);  \
                             \
    element_mul(t0, b, Zy);  \
    element_mul(c, a, Zx);   \
    element_add(c, c, t0);   \
    element_neg(c, c);       \
                             \
    f_miller_evalfn();       \
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
                            \
    f_miller_evalfn();      \
  }

  element_init(a, Px->field);
  element_init(b, a->field);
  element_init(c, a->field);
  element_init(t0, a->field);
  element_init(e0, res->field);
  element_init(e1, Qx->field);

  element_init(v, res->field);
  element_init(Z, P->field);

  element_set(Z, P);
  Zx = curve_x_coord(Z);
  Zy = curve_y_coord(Z);

  element_set1(v);
  m = mpz_sizeinbase(q, 2) - 2;

  //TODO: sliding NAF
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
  element_clear(e1);
  #undef do_term
  #undef f_miller_evalfn
  #undef do_tangent
  #undef do_line
}

static void f_tateexp(element_t out) {
  element_t x, y, epow;
  f_pairing_data_ptr p = out->field->pairing->data;
  element_init(x, p->Fq12);
  element_init(y, p->Fq12);
  element_init(epow, p->Fq2);

  #define qpower(e1, e) {                                         \
    element_set(element_item(e1, 0), element_item(out, 0));       \
    element_mul(element_item(e1, 1), element_item(out, 1), e);    \
    element_square(epow, e);                                      \
    element_mul(element_item(e1, 2), element_item(out, 2), epow); \
    element_mul(epow, epow, e);                                   \
    element_mul(element_item(e1, 3), element_item(out, 3), epow); \
    element_mul(epow, epow, e);                                   \
    element_mul(element_item(e1, 4), element_item(out, 4), epow); \
    element_mul(epow, epow, e);                                   \
    element_mul(element_item(e1, 5), element_item(out, 5), epow); \
  }

  qpower(y, p->xpowq8);
  qpower(x, p->xpowq6);
  element_mul(y, y, x);
  qpower(x, p->xpowq2);
  element_mul(x, x, out);
  element_invert(x, x);
  element_mul(out, y, x);

  element_clear(epow);
  element_clear(x);
  element_clear(y);
  element_pow_mpz(out, out, p->tateexp);
  #undef qpower
}

static void f_finalpow(element_t out) {
  f_tateexp(out->data);
}

static void f_pairing(element_ptr out, element_ptr in1, element_ptr in2,
    pairing_t pairing) {
  element_ptr Qbase = in2;
  element_t x, y;
  f_pairing_data_ptr p = pairing->data;

  element_init(x, p->Fq2);
  element_init(y, p->Fq2);
  //map from twist: (x, y) --> (v^-2 x, v^-3 y)
  //where v is the sixth root used to construct the twist
  //i.e. v^6 = -alpha
  //thus v^-2 = -alpha^-1 v^4
  //and  v^-3 = -alpha^-1 v^3
  element_mul(x, curve_x_coord(Qbase), p->negalphainv);
  element_mul(y, curve_y_coord(Qbase), p->negalphainv);

  cc_miller_no_denom(out, pairing->r, in1, x, y, p->negalpha);

  element_clear(x);
  element_clear(y);

  f_tateexp(out);
}

static void f_pairing_clear(pairing_t pairing) {
  field_clear(pairing->GT);
  f_pairing_data_ptr p = pairing->data;
  element_clear(p->negalpha);
  element_clear(p->negalphainv);
  mpz_clear(p->tateexp);
  element_clear(p->xpowq2);
  element_clear(p->xpowq6);
  element_clear(p->xpowq8);
  field_clear(p->Etwist);
  field_clear(p->Eq);

  field_clear(p->Fq12);
  field_clear(p->Fq2x);
  field_clear(p->Fq2);
  field_clear(p->Fq);
  pbc_free(p);

  mpz_clear(pairing->r);
  field_clear(pairing->Zr);
}

static void f_init_pairing(pairing_t pairing, void *data) {
  f_param_ptr param = data;
  f_pairing_data_ptr p;
  element_t irred;
  element_t e0, e1, e2;
  p = pairing->data = pbc_malloc(sizeof(f_pairing_data_t));
  mpz_init(pairing->r);
  mpz_set(pairing->r, param->r);
  field_init_fp(pairing->Zr, pairing->r);
  field_init_fp(p->Fq, param->q);
  p->Fq->nqr = pbc_malloc(sizeof(element_t));
  element_init(p->Fq->nqr, p->Fq);
  element_set_mpz(p->Fq->nqr, param->beta);
  field_init_quadratic(p->Fq2, p->Fq);
  field_init_poly(p->Fq2x, p->Fq2);
  element_init(irred, p->Fq2x);
  // Call poly_set_coeff1() first so we can use element_item() for the other
  // coefficients.
  poly_set_coeff1(irred, 6);

  element_init(p->negalpha, p->Fq2);
  element_init(p->negalphainv, p->Fq2);
  element_set_mpz(element_x(p->negalpha), param->alpha0);
  element_set_mpz(element_y(p->negalpha), param->alpha1);

  element_set(element_item(irred, 0), p->negalpha);
  field_init_polymod(p->Fq12, irred);
  element_neg(p->negalpha, p->negalpha);
  element_invert(p->negalphainv, p->negalpha);
  element_clear(irred);

  element_init(e0, p->Fq);
  element_init(e1, p->Fq);
  element_init(e2, p->Fq2);

  // Initialize the curve Y^2 = X^3 + b.
  element_set_mpz(e1, param->b);
  field_init_curve_ab(p->Eq, e0, e1, pairing->r, NULL);

  // Initialize the curve Y^2 = X^3 - alpha0 b - alpha1 sqrt(beta) b.
  element_set_mpz(e0, param->alpha0);
  element_neg(e0, e0);
  element_mul(element_x(e2), e0, e1);
  element_set_mpz(e0, param->alpha1);
  element_neg(e0, e0);
  element_mul(element_y(e2), e0, e1);
  element_clear(e0);
  element_init(e0, p->Fq2);
  field_init_curve_ab(p->Etwist, e0, e2, pairing->r, NULL);
  element_clear(e0);
  element_clear(e1);
  element_clear(e2);

  mpz_t ndonr;
  mpz_init(ndonr);
  // ndonr temporarily holds the trace.
  mpz_sub(ndonr, param->q, param->r);
  mpz_add_ui(ndonr, ndonr, 1);
  // TODO: We can use a smaller quotient_cmp, but I have to figure out
  // BN curves again.
  pbc_mpz_curve_order_extn(ndonr, param->q, ndonr, 12);
  mpz_divexact(ndonr, ndonr, param->r);
  mpz_divexact(ndonr, ndonr, param->r);
  field_curve_set_quotient_cmp(p->Etwist, ndonr);
  mpz_clear(ndonr);

  pairing->G1 = p->Eq;
  pairing->G2 = p->Etwist;
  pairing_GT_init(pairing, p->Fq12);
  pairing->finalpow = f_finalpow;
  pairing->map = f_pairing;
  pairing->clear_func = f_pairing_clear;

  mpz_init(p->tateexp);
  /* unoptimized tate exponent
  mpz_pow_ui(p->tateexp, param->q, 12);
  mpz_sub_ui(p->tateexp, p->tateexp, 1);
  mpz_divexact(p->tateexp, p->tateexp, param->r);
  */
  mpz_ptr z = p->tateexp;
  mpz_mul(z, param->q, param->q);
  mpz_sub_ui(z, z, 1);
  mpz_mul(z, z, param->q);
  mpz_mul(z, z, param->q);
  mpz_add_ui(z, z, 1);
  mpz_divexact(z, z, param->r);

  element_init(p->xpowq2, p->Fq2);
  element_init(p->xpowq6, p->Fq2);
  element_init(p->xpowq8, p->Fq2);
  element_t xpowq;
  element_init(xpowq, p->Fq12);

  //there are smarter ways since we know q = 1 mod 6
  //and that x^6 = -alpha
  //but this is fast enough
  element_set1(element_item(xpowq, 1));
  element_pow_mpz(xpowq, xpowq, param->q);
  element_pow_mpz(xpowq, xpowq, param->q);
  element_set(p->xpowq2, element_item(xpowq, 1));

  element_pow_mpz(xpowq, xpowq, param->q);
  element_pow_mpz(xpowq, xpowq, param->q);
  element_pow_mpz(xpowq, xpowq, param->q);
  element_pow_mpz(xpowq, xpowq, param->q);
  element_set(p->xpowq6, element_item(xpowq, 1));

  element_pow_mpz(xpowq, xpowq, param->q);
  element_pow_mpz(xpowq, xpowq, param->q);
  element_set(p->xpowq8, element_item(xpowq, 1));

  element_clear(xpowq);
}

static void f_init(pbc_param_ptr p) {
  static pbc_param_interface_t interface = {{
    f_clear,
    f_init_pairing,
    f_out_str,
  }};
  p->api = interface;
  f_param_ptr fp = p->data = pbc_malloc(sizeof(*fp));
  mpz_init(fp->q);
  mpz_init(fp->r);
  mpz_init(fp->b);
  mpz_init(fp->beta);
  mpz_init(fp->alpha0);
  mpz_init(fp->alpha1);
}

// Public interface:

int pbc_param_init_f(pbc_param_ptr par, struct symtab_s *tab) {
  f_init(par);
  f_param_ptr p = par->data;

  int err = 0;
  err += lookup_mpz(p->q, tab, "q");
  err += lookup_mpz(p->r, tab, "r");
  err += lookup_mpz(p->b, tab, "b");
  err += lookup_mpz(p->beta, tab, "beta");
  err += lookup_mpz(p->alpha0, tab, "alpha0");
  err += lookup_mpz(p->alpha1, tab, "alpha1");
  return err;
}

void pbc_param_init_f_gen(pbc_param_t p, int bits) {
  f_init(p);
  f_param_ptr fp = p->data;
  //36 is a 6-bit number
  int xbit = (bits - 6) / 4;
  //TODO: use binary search to find smallest appropriate x
  mpz_t x, t;
  mpz_ptr q = fp->q;
  mpz_ptr r = fp->r;
  mpz_ptr b = fp->b;
  field_t Fq, Fq2, Fq2x;
  element_t e1;
  element_t f;
  field_t c;
  element_t P;

  mpz_init(x);
  mpz_init(t);
  mpz_setbit(x, xbit);
  for (;;) {
    mpz_mul(t, x, x);
    mpz_mul_ui(t, t, 6);
    mpz_add_ui(t, t, 1);
    tryminusx(q, x);
    mpz_sub(r, q, t);
    mpz_add_ui(r, r, 1);
    if (mpz_probab_prime_p(q, 10) && mpz_probab_prime_p(r, 10)) break;

    tryplusx(q, x);
    mpz_sub(r, q, t);
    mpz_add_ui(r, r, 1);
    if (mpz_probab_prime_p(q, 10) && mpz_probab_prime_p(r, 10)) break;

    mpz_add_ui(x, x, 1);
  }

  field_init_fp(Fq, q);
  element_init(e1, Fq);

  for (;;) {
    element_random(e1);
    field_init_curve_b(c, e1, r, NULL);
    element_init(P, c);

    element_random(P);

    element_mul_mpz(P, P, r);
    if (element_is0(P)) break;
    element_clear(P);
    field_clear(c);
  }
  element_to_mpz(b, e1);
  element_clear(e1);
  field_init_quadratic(Fq2, Fq);
  element_to_mpz(fp->beta, field_get_nqr(Fq));
  field_init_poly(Fq2x, Fq2);
  element_init(f, Fq2x);

  // Find an irreducible polynomial of the form f = x^6 + alpha.
  // Call poly_set_coeff1() first so we can use element_item() for the other
  // coefficients.
  poly_set_coeff1(f, 6);
  for (;;) {
    element_random(element_item(f, 0));
    if (poly_is_irred(f)) break;
  }

  //extend F_q^2 using f = x^6 + alpha
  //see if sextic twist contains a subgroup of order r
  //if not, it's the wrong twist: replace alpha with alpha^5
  {
    field_t ctest;
    element_t Ptest;
    mpz_t z0, z1;
    mpz_init(z0);
    mpz_init(z1);
    element_init(e1, Fq2);
    element_set_mpz(e1, fp->b);
    element_mul(e1, e1, element_item(f, 0));
    element_neg(e1, e1);

    field_init_curve_b(ctest, e1, r, NULL);
    element_init(Ptest, ctest);
    element_random(Ptest);

    //I'm not sure what the #E'(F_q^2) is, but
    //it definitely divides n_12 = #E(F_q^12). It contains a
    //subgroup of order r if and only if
    //(n_12 / r^2)P != O for some (in fact most) P in E'(F_q^6)
    mpz_pow_ui(z0, q, 12);
    mpz_add_ui(z0, z0, 1);
    pbc_mpz_trace_n(z1, q, t, 12);
    mpz_sub(z1, z0, z1);
    mpz_mul(z0, r, r);
    mpz_divexact(z1, z1, z0);

    element_mul_mpz(Ptest, Ptest, z1);
    if (element_is0(Ptest)) {
      mpz_set_ui(z0, 5);
      element_pow_mpz(element_item(f, 0), element_item(f, 0), z0);
    }
    element_clear(e1);
    element_clear(Ptest);
    field_clear(ctest);
    mpz_clear(z0);
    mpz_clear(z1);
  }

  element_to_mpz(fp->alpha0, element_x(element_item(f, 0)));
  element_to_mpz(fp->alpha1, element_y(element_item(f, 0)));

  element_clear(f);

  field_clear(Fq2x);
  field_clear(Fq2);
  field_clear(Fq);

  mpz_clear(t);
  mpz_clear(x);
}
