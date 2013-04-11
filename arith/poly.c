#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdint.h> // for intptr_t
#include <stdlib.h>
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_field.h"
#include "pbc_multiz.h"
#include "pbc_poly.h"
#include "pbc_memory.h"
#include "misc/darray.h"

// == Polynomial rings ==
//
// Per-field data:
typedef struct {
  field_ptr field;  // Ring where coefficients live.
  fieldmap mapbase; // Map element from underlying field to constant term.
} *pfptr;

// Per-element data:
//TODO: Would we ever need any field besides coeff?
typedef struct {
  // The coefficients are held in a darray which is resized as needed.
  // The last array entry represents the leading coefficient and should be
  // nonzero. An empty darray represents 0.
  darray_t coeff;
} *peptr;

// == Polynomial modulo rings ==
//
// Per-field data:
typedef struct {
  field_ptr field;   // Base field.
  fieldmap mapbase;  // Similar to mapbase above.
  int n;             // Degree of extension.
  element_t poly;    // Polynomial of degree n.
  element_t *xpwr;   // x^n,...,x^{2n-2} mod poly
} *mfptr;
// Per-element data: just a pointer to an array of element_t. This array always
// has size n.

// Add or remove coefficients until there are exactly n of them. Any new
// coefficients are initialized to zero, which violates the invariant that the
// leading coefficient must be nonzero. Thus routines calling this function
// must check for this and fix the polynomial if necessary, e.g. by calling
// poly_remove_leading_zeroes().
static void poly_alloc(element_ptr e, int n) {
  pfptr pdp = e->field->data;
  peptr p = e->data;
  element_ptr e0;
  int k = p->coeff->count;
  while (k < n) {
    e0 = pbc_malloc(sizeof(element_t));
    element_init(e0, pdp->field);
    darray_append(p->coeff, e0);
    k++;
  }
  while (k > n) {
    k--;
    e0 = darray_at(p->coeff, k);
    element_clear(e0);
    pbc_free(e0);
    darray_remove_last(p->coeff);
  }
}

static void poly_init(element_ptr e) {
  peptr p = e->data = pbc_malloc(sizeof(*p));
  darray_init(p->coeff);
}

static void poly_clear(element_ptr e) {
  peptr p = e->data;

  poly_alloc(e, 0);
  darray_clear(p->coeff);
  pbc_free(e->data);
}

// Some operations may zero a leading coefficient, which will cause other
// routines to fail. After such an operation, this function should be called,
// as it strips all leading zero coefficients and frees the memory they
// occupied, reestablishing the guarantee that the last element of the array
// is nonzero.
static void poly_remove_leading_zeroes(element_ptr e) {
  peptr p = e->data;
  int n = p->coeff->count - 1;
  while (n >= 0) {
    element_ptr e0 = p->coeff->item[n];
    if (!element_is0(e0)) return;
    element_clear(e0);
    pbc_free(e0);
    darray_remove_last(p->coeff);
    n--;
  }
}

static void poly_set0(element_ptr e) {
  poly_alloc(e, 0);
}

static void poly_set1(element_ptr e) {
  peptr p = e->data;
  element_ptr e0;

  poly_alloc(e, 1);
  e0 = p->coeff->item[0];
  element_set1(e0);
}

static int poly_is0(element_ptr e) {
  peptr p = e->data;
  return !p->coeff->count;
}

static int poly_is1(element_ptr e) {
  peptr p = e->data;
  if (p->coeff->count == 1) {
    return element_is1(p->coeff->item[0]);
  }
  return 0;
}

static void poly_set_si(element_ptr e, signed long int op) {
  peptr p = e->data;
  element_ptr e0;

  poly_alloc(e, 1);
  e0 = p->coeff->item[0];
  element_set_si(e0, op);
  poly_remove_leading_zeroes(e);
}

static void poly_set_mpz(element_ptr e, mpz_ptr op) {
  peptr p = e->data;

  poly_alloc(e, 1);
  element_set_mpz(p->coeff->item[0], op);
  poly_remove_leading_zeroes(e);
}

static void poly_set_multiz(element_ptr e, multiz op) {
  if (multiz_is_z(op)) {
    // TODO: Remove unnecessary copy.
    mpz_t z;
    mpz_init(z);
    multiz_to_mpz(z, op);
    poly_set_mpz(e, z);
    mpz_clear(z);
    return;
  }
  peptr p = e->data;
  int n = multiz_count(op);
  poly_alloc(e, n);
  int i;
  for(i = 0; i < n; i++) {
    element_set_multiz(p->coeff->item[i], multiz_at(op, i));
  }
  poly_remove_leading_zeroes(e);
}

static void poly_set(element_ptr dst, element_ptr src) {
  peptr psrc = src->data;
  peptr pdst = dst->data;
  int i;

  poly_alloc(dst, psrc->coeff->count);
  for (i=0; i<psrc->coeff->count; i++) {
    element_set(pdst->coeff->item[i], psrc->coeff->item[i]);
  }
}

static int poly_coeff_count(element_ptr e) {
  return ((peptr) e->data)->coeff->count;
}

static element_ptr poly_coeff(element_ptr e, int n) {
  peptr ep = e->data;
  PBC_ASSERT(n < poly_coeff_count(e), "coefficient out of range");
  return (element_ptr) ep->coeff->item[n];
}

static int poly_sgn(element_ptr f) {
  int res = 0;
  int i;
  int n = poly_coeff_count(f);
  for (i=0; i<n; i++) {
    res = element_sgn(poly_coeff(f, i));
    if (res) break;
  }
  return res;
}

static void poly_add(element_ptr sum, element_ptr f, element_ptr g) {
  int i, n, n1;
  element_ptr big;

  n = poly_coeff_count(f);
  n1 = poly_coeff_count(g);
  if (n > n1) {
    big = f;
    n = n1;
    n1 = poly_coeff_count(f);
  } else {
    big = g;
  }

  poly_alloc(sum, n1);
  for (i=0; i<n; i++) {
    element_add(poly_coeff(sum, i), poly_coeff(f, i), poly_coeff(g, i));
  }
  for (; i<n1; i++) {
    element_set(poly_coeff(sum, i), poly_coeff(big, i));
  }
  poly_remove_leading_zeroes(sum);
}

static void poly_sub(element_ptr diff, element_ptr f, element_ptr g) {
  int i, n, n1;
  element_ptr big;

  n = poly_coeff_count(f);
  n1 = poly_coeff_count(g);
  if (n > n1) {
    big = f;
    n = n1;
    n1 = poly_coeff_count(f);
  } else {
    big = g;
  }

  poly_alloc(diff, n1);
  for (i=0; i<n; i++) {
    element_sub(poly_coeff(diff, i), poly_coeff(f, i), poly_coeff(g, i));
  }
  for (; i<n1; i++) {
    if (big == f) {
      element_set(poly_coeff(diff, i), poly_coeff(big, i));
    } else {
      element_neg(poly_coeff(diff, i), poly_coeff(big, i));
    }
  }
  poly_remove_leading_zeroes(diff);
}

static void poly_neg(element_ptr f, element_ptr g) {
  peptr pf = f->data;
  peptr pg = g->data;
  int i, n;

  n = pg->coeff->count;
  poly_alloc(f, n);
  for (i=0; i<n; i++) {
    element_neg(pf->coeff->item[i], pg->coeff->item[i]);
  }
}

static void poly_double(element_ptr f, element_ptr g) {
  peptr pf = f->data;
  peptr pg = g->data;
  int i, n;

  n = pg->coeff->count;
  poly_alloc(f, n);
  for (i=0; i<n; i++) {
    element_double(pf->coeff->item[i], pg->coeff->item[i]);
  }
}

static void poly_mul_mpz(element_ptr f, element_ptr g, mpz_ptr z) {
  peptr pf = f->data;
  peptr pg = g->data;
  int i, n;

  n = pg->coeff->count;
  poly_alloc(f, n);
  for (i=0; i<n; i++) {
    element_mul_mpz(pf->coeff->item[i], pg->coeff->item[i], z);
  }
}

static void poly_mul_si(element_ptr f, element_ptr g, signed long int z) {
  peptr pf = f->data;
  peptr pg = g->data;
  int i, n;

  n = pg->coeff->count;
  poly_alloc(f, n);
  for (i=0; i<n; i++) {
    element_mul_si(pf->coeff->item[i], pg->coeff->item[i], z);
  }
}

static void poly_mul(element_ptr r, element_ptr f, element_ptr g) {
  peptr pprod;
  peptr pf = f->data;
  peptr pg = g->data;
  pfptr pdp = r->field->data;
  int fcount = pf->coeff->count;
  int gcount = pg->coeff->count;
  int i, j, n;
  element_t prod;
  element_t e0;

  if (!fcount || !gcount) {
    element_set0(r);
    return;
  }
  element_init(prod, r->field);
  pprod = prod->data;
  n = fcount + gcount - 1;
  poly_alloc(prod, n);
  element_init(e0, pdp->field);
  for (i=0; i<n; i++) {
    element_ptr x = pprod->coeff->item[i];
    element_set0(x);
    for (j=0; j<=i; j++) {
      if (j < fcount && i - j < gcount) {
        element_mul(e0, pf->coeff->item[j], pg->coeff->item[i - j]);
        element_add(x, x, e0);
      }
    }
  }
  poly_remove_leading_zeroes(prod);
  element_set(r, prod);
  element_clear(e0);
  element_clear(prod);
}

static void polymod_random(element_ptr e) {
  element_t *coeff = e->data;
  int i, n = polymod_field_degree(e->field);

  for (i=0; i<n; i++) {
    element_random(coeff[i]);
  }
}

static void polymod_from_hash(element_ptr e, void *data, int len) {
  // TODO: Improve this.
  element_t *coeff = e->data;
  int i, n = polymod_field_degree(e->field);
  for (i=0; i<n; i++) {
    element_from_hash(coeff[i], data, len);
  }
}

static size_t poly_out_str(FILE *stream, int base, element_ptr e) {
  int i;
  int n = poly_coeff_count(e);
  size_t result = 2, status;

  /*
  if (!n) {
    if (EOF == fputs("[0]", stream)) return 0;
    return 3;
  }
  */
  if (EOF == fputc('[', stream)) return 0;
  for (i=0; i<n; i++) {
    if (i) {
      if (EOF == fputs(", ", stream)) return 0;
      result += 2;
    }
    status = element_out_str(stream, base, poly_coeff(e, i));
    if (!status) return 0;
    result += status;
  }
  if (EOF == fputc(']', stream)) return 0;
  return result;
}

static int poly_snprint(char *s, size_t size, element_ptr e) {
  int i;
  int n = poly_coeff_count(e);
  size_t result = 0, left;
  int status;

  #define clip_sub() {                         \
    result += status;                          \
    left = result >= size ? 0 : size - result; \
  }

  status = snprintf(s, size, "[");
  if (status < 0) return status;
  clip_sub();

  for (i=0; i<n; i++) {
    if (i) {
      status = snprintf(s + result, left, ", ");
      if (status < 0) return status;
      clip_sub();
    }
    status = element_snprint(s + result, left, poly_coeff(e, i));
    if (status < 0) return status;
    clip_sub();
  }
  status = snprintf(s + result, left, "]");
  if (status < 0) return status;
  return result + status;
  #undef clip_sub
}

static void poly_div(element_ptr quot, element_ptr rem,
    element_ptr a, element_ptr b) {
  peptr pq, pr;
  pfptr pdp = a->field->data;
  element_t q, r;
  element_t binv, e0;
  element_ptr qe;
  int m, n;
  int i, k;

  if (element_is0(b)) pbc_die("division by zero");
  n = poly_degree(b);
  m = poly_degree(a);
  if (n > m) {
    element_set(rem, a);
    element_set0(quot);
    return;
  }
  element_init(r, a->field);
  element_init(q, a->field);
  element_init(binv, pdp->field);
  element_init(e0, pdp->field);
  pq = q->data;
  pr = r->data;
  element_set(r, a);
  k = m - n;
  poly_alloc(q, k + 1);
  element_invert(binv, poly_coeff(b, n));
  while (k >= 0) {
    qe = pq->coeff->item[k];
    element_mul(qe, binv, pr->coeff->item[m]);
    for (i=0; i<=n; i++) {
      element_mul(e0, qe, poly_coeff(b, i));
      element_sub(pr->coeff->item[i + k], pr->coeff->item[i + k], e0);
    }
    k--;
    m--;
  }
  poly_remove_leading_zeroes(r);
  element_set(quot, q);
  element_set(rem, r);

  element_clear(q);
  element_clear(r);
  element_clear(e0);
  element_clear(binv);
}

static void poly_invert(element_ptr res, element_ptr f, element_ptr m) {
  element_t q, r0, r1, r2;
  element_t b0, b1, b2;
  element_t inv;

  element_init(b0, res->field);
  element_init(b1, res->field);
  element_init(b2, res->field);
  element_init(q, res->field);
  element_init(r0, res->field);
  element_init(r1, res->field);
  element_init(r2, res->field);
  element_init(inv, poly_base_field(res));
  element_set0(b0);
  element_set1(b1);
  element_set(r0, m);
  element_set(r1, f);

  for (;;) {
    poly_div(q, r2, r0, r1);
    if (element_is0(r2)) break;
    element_mul(b2, b1, q);
    element_sub(b2, b0, b2);
    element_set(b0, b1);
    element_set(b1, b2);
    element_set(r0, r1);
    element_set(r1, r2);
  }
  element_invert(inv, poly_coeff(r1, 0));
  poly_const_mul(res, inv, b1);
  element_clear(inv);
  element_clear(q);
  element_clear(r0);
  element_clear(r1);
  element_clear(r2);
  element_clear(b0);
  element_clear(b1);
  element_clear(b2);
}

static void poly_to_polymod_truncate(element_ptr e, element_ptr f) {
  mfptr p = e->field->data;
  element_t *coeff = e->data;
  int i;
  int n;
  n = poly_coeff_count(f);
  if (n > p->n) n = p->n;

  for (i=0; i<n; i++) {
    element_set(coeff[i], poly_coeff(f, i));
  }
  for (; i<p->n; i++) {
    element_set0(coeff[i]);
  }
}

static void polymod_to_poly(element_ptr f, element_ptr e) {
  mfptr p = e->field->data;
  element_t *coeff = e->data;
  int i, n = p->n;
  poly_alloc(f, n);
  for (i=0; i<n; i++) {
    element_set(poly_coeff(f, i), coeff[i]);
  }
  poly_remove_leading_zeroes(f);
}

static void polymod_invert(element_ptr r, element_ptr e) {
  mfptr p = r->field->data;
  element_ptr minpoly = p->poly;
  element_t f, r1;

  element_init(f, minpoly->field);
  element_init(r1, minpoly->field);
  polymod_to_poly(f, e);

  poly_invert(r1, f, p->poly);

  poly_to_polymod_truncate(r, r1);

  element_clear(f);
  element_clear(r1);
}

static int poly_cmp(element_ptr f, element_ptr g) {
  int i;
  int n = poly_coeff_count(f);
  int n1 = poly_coeff_count(g);
  if (n != n1) return 1;
  for (i=0; i<n; i++) {
    if (element_cmp(poly_coeff(f, i), poly_coeff(g, i))) return 1;
  }
  return 0;
}

static void field_clear_poly(field_ptr f) {
  pfptr p = f->data;
  pbc_free(p);
}

// 2 bytes hold the number of terms, then the terms follow.
// Bad for sparse polynomials.
static int poly_length_in_bytes(element_t p) {
  int count = poly_coeff_count(p);
  int result = 2;
  int i;
  for (i=0; i<count; i++) {
    result += element_length_in_bytes(poly_coeff(p, i));
  }
  return result;
}

static int poly_to_bytes(unsigned char *buf, element_t p) {
  int count = poly_coeff_count(p);
  int result = 2;
  int i;
  buf[0] = (unsigned char) count;
  buf[1] = (unsigned char) (count >> 8);
  for (i=0; i<count; i++) {
    result += element_to_bytes(&buf[result], poly_coeff(p, i));
  }
  return result;
}

static int poly_from_bytes(element_t p, unsigned char *buf) {
  int result = 2;
  int count = buf[0] + buf[1] * 256;
  int i;
  poly_alloc(p, count);
  for (i=0; i<count; i++) {
    result += element_from_bytes(poly_coeff(p, i), &buf[result]);
  }
  return result;
}

// Is this useful? This returns to_mpz(constant term).
static void poly_to_mpz(mpz_t z, element_ptr e) {
  if (!poly_coeff_count(e)) {
    mpz_set_ui(z, 0);
  } else {
    element_to_mpz(z, poly_coeff(e, 0));
  }
}

static void poly_out_info(FILE *str, field_ptr f) {
  pfptr p = f->data;
  fprintf(str, "Polynomial ring over ");
  field_out_info(str, p->field);
}

static void field_clear_polymod(field_ptr f) {
  mfptr p = f->data;
  int i, n = p->n;

  for (i=0; i<n; i++) {
    element_clear(p->xpwr[i]);
  }
  pbc_free(p->xpwr);

  element_clear(p->poly);
  pbc_free(f->data);
}

static int polymod_is_sqr(element_ptr e) {
  int res;
  mpz_t z;
  element_t e0;

  element_init(e0, e->field);
  mpz_init(z);
  mpz_sub_ui(z, e->field->order, 1);
  mpz_divexact_ui(z, z, 2);

  element_pow_mpz(e0, e, z);
  res = element_is1(e0);
  element_clear(e0);
  mpz_clear(z);
  return res;
}

// Find a square root in a polynomial modulo ring using Cantor-Zassenhaus aka
// Legendre's method.
static void polymod_sqrt(element_ptr res, element_ptr a) {
  // TODO: Use a faster method? See Bernstein.
  field_t kx;
  element_t f;
  element_t r, s;
  element_t e0;
  mpz_t z;

  field_init_poly(kx, a->field);
  mpz_init(z);
  element_init(f, kx);
  element_init(r, kx);
  element_init(s, kx);
  element_init(e0, a->field);

  poly_alloc(f, 3);
  element_set1(poly_coeff(f, 2));
  element_neg(poly_coeff(f, 0), a);

  mpz_sub_ui(z, a->field->order, 1);
  mpz_divexact_ui(z, z, 2);
  for (;;) {
    int i;
    element_ptr x;
    element_ptr e1, e2;

    poly_alloc(r, 2);
    element_set1(poly_coeff(r, 1));
    x = poly_coeff(r, 0);
    element_random(x);
    element_mul(e0, x, x);
    if (!element_cmp(e0, a)) {
      element_set(res, x);
      break;
    }
    element_set1(s);
    //TODO: this can be optimized greatly
    //since we know r has the form ax + b
    for (i = mpz_sizeinbase(z, 2) - 1; i >=0; i--) {
      element_mul(s, s, s);
      if (poly_degree(s) == 2) {
        e1 = poly_coeff(s, 0);
        e2 = poly_coeff(s, 2);
        element_mul(e0, e2, a);
        element_add(e1, e1, e0);
        poly_alloc(s, 2);
        poly_remove_leading_zeroes(s);
      }
      if (mpz_tstbit(z, i)) {
        element_mul(s, s, r);
        if (poly_degree(s) == 2) {
          e1 = poly_coeff(s, 0);
          e2 = poly_coeff(s, 2);
          element_mul(e0, e2, a);
          element_add(e1, e1, e0);
          poly_alloc(s, 2);
          poly_remove_leading_zeroes(s);
        }
      }
    }
    if (poly_degree(s) < 1) continue;
    element_set1(e0);
    e1 = poly_coeff(s, 0);
    e2 = poly_coeff(s, 1);
    element_add(e1, e1, e0);
    element_invert(e0, e2);
    element_mul(e0, e0, e1);
    element_mul(e2, e0, e0);
    if (!element_cmp(e2, a)) {
      element_set(res, e0);
      break;
    }
  }

  mpz_clear(z);
  element_clear(f);
  element_clear(r);
  element_clear(s);
  element_clear(e0);
  field_clear(kx);
}

static int polymod_to_bytes(unsigned char *data, element_t f) {
  mfptr p = f->field->data;
  element_t *coeff = f->data;
  int i, n = p->n;
  int len = 0;
  for (i=0; i<n; i++) {
    len += element_to_bytes(data + len, coeff[i]);
  }
  return len;
}

static int polymod_length_in_bytes(element_t f) {
  mfptr p = f->field->data;
  element_t *coeff = f->data;
  int i, n = p->n;
  int res = 0;

  for (i=0; i<n; i++) {
    res += element_length_in_bytes(coeff[i]);
  }

  return res;
}

static int polymod_from_bytes(element_t f, unsigned char *data) {
  mfptr p = f->field->data;
  element_t *coeff = f->data;
  int i, n = p->n;
  int len = 0;

  for (i=0; i<n; i++) {
    len += element_from_bytes(coeff[i], data + len);
  }
  return len;
}

static void polymod_init(element_t e) {
  int i;
  mfptr p = e->field->data;
  int n = p->n;
  element_t *coeff;
  coeff = e->data = pbc_malloc(sizeof(element_t) * n);

  for (i=0; i<n; i++) {
    element_init(coeff[i], p->field);
  }
}

static void polymod_clear(element_t e) {
  mfptr p = e->field->data;
  element_t *coeff = e->data;
  int i, n = p->n;
  for (i=0; i<n; i++) {
    element_clear(coeff[i]);
  }
  pbc_free(e->data);
}

static void polymod_set_si(element_t e, signed long int x) {
  mfptr p = e->field->data;
  element_t *coeff = e->data;
  int i, n = p->n;
  element_set_si(coeff[0], x);
  for (i=1; i<n; i++) {
    element_set0(coeff[i]);
  }
}

static void polymod_set_mpz(element_t e, mpz_t z) {
  mfptr p = e->field->data;
  element_t *coeff = e->data;
  int i, n = p->n;
  element_set_mpz(coeff[0], z);
  for (i=1; i<n; i++) {
    element_set0(coeff[i]);
  }
}

static void polymod_set(element_t e, element_t f) {
  mfptr p = e->field->data;
  element_t *dst = e->data, *src = f->data;
  int i, n = p->n;
  for (i=0; i<n; i++) {
    element_set(dst[i], src[i]);
  }
}

static void polymod_neg(element_t e, element_t f) {
  mfptr p = e->field->data;
  element_t *dst = e->data, *src = f->data;
  int i, n = p->n;
  for (i=0; i<n; i++) {
    element_neg(dst[i], src[i]);
  }
}

static int polymod_cmp(element_ptr f, element_ptr g) {
  mfptr p = f->field->data;
  element_t *c1 = f->data, *c2 = g->data;
  int i, n = p->n;
  for (i=0; i<n; i++) {
    if (element_cmp(c1[i], c2[i])) return 1;
  }
  return 0;
}

static void polymod_add(element_t r, element_t e, element_t f) {
  mfptr p = r->field->data;
  element_t *dst = r->data, *s1 = e->data, *s2 = f->data;
  int i, n = p->n;
  for (i=0; i<n; i++) {
    element_add(dst[i], s1[i], s2[i]);
  }
}

static void polymod_double(element_t r, element_t f) {
  mfptr p = r->field->data;
  element_t *dst = r->data, *s1 = f->data;
  int i, n = p->n;
  for (i=0; i<n; i++) {
    element_double(dst[i], s1[i]);
  }
}

static void polymod_sub(element_t r, element_t e, element_t f) {
  mfptr p = r->field->data;
  element_t *dst = r->data, *s1 = e->data, *s2 = f->data;
  int i, n = p->n;
  for (i=0; i<n; i++) {
    element_sub(dst[i], s1[i], s2[i]);
  }
}

static void polymod_mul_mpz(element_t e, element_t f, mpz_ptr z) {
  mfptr p = e->field->data;
  element_t *dst = e->data, *src = f->data;
  int i, n = p->n;
  for (i=0; i<n; i++) {
    element_mul_mpz(dst[i], src[i], z);
  }
}

static void polymod_mul_si(element_t e, element_t f, signed long int z) {
  mfptr p = e->field->data;
  element_t *dst = e->data, *src = f->data;
  int i, n = p->n;
  for (i=0; i<n; i++) {
    element_mul_si(dst[i], src[i], z);
  }
}

// Karatsuba multiplication for degree 2 polynomials.
static void kar_poly_2(element_t *dst, element_t c3, element_t c4, element_t *s1, element_t *s2, element_t *scratch) {
  element_ptr c01, c02, c12;

  c12 = scratch[0];
  c02 = scratch[1];
  c01 = scratch[2];

  element_add(c3, s1[0], s1[1]);
  element_add(c4, s2[0], s2[1]);
  element_mul(c01, c3, c4);
  element_add(c3, s1[0], s1[2]);
  element_add(c4, s2[0], s2[2]);
  element_mul(c02, c3, c4);
  element_add(c3, s1[1], s1[2]);
  element_add(c4, s2[1], s2[2]);
  element_mul(c12, c3, c4);

  element_mul(dst[1], s1[1], s2[1]);

  // Constant term.
  element_mul(dst[0], s1[0], s2[0]);

  // Coefficient of x^4.
  element_mul(c4, s1[2], s2[2]);

  // Coefficient of x^3.
  element_add(c3, dst[1], c4);
  element_sub(c3, c12, c3);

  // Coefficient of x^2.
  element_add(dst[2], c4, dst[0]);
  element_sub(c02, c02, dst[2]);
  element_add(dst[2], dst[1], c02);

  // Coefficient of x.
  element_sub(c01, c01, dst[0]);
  element_sub(dst[1], c01, dst[1]);
}

// Degree 3, 6 polynomial moduli have dedicated routines for multiplication.
static void polymod_mul_degree3(element_ptr res, element_ptr e, element_ptr f) {
  mfptr p = res->field->data;
  element_t *dst = res->data, *s1 = e->data, *s2 = f->data;
  element_t c3, c4;
  element_t p0;

  element_init(p0, res->field);
  element_init(c3, p->field);
  element_init(c4, p->field);

  kar_poly_2(dst, c3, c4, s1, s2, p0->data);

  polymod_const_mul(p0, c3, p->xpwr[0]);
  element_add(res, res, p0);
  polymod_const_mul(p0, c4, p->xpwr[1]);
  element_add(res, res, p0);

  element_clear(p0);
  element_clear(c3);
  element_clear(c4);
}

static void polymod_mul_degree6(element_ptr res, element_ptr e, element_ptr f) {
  mfptr p = res->field->data;
  element_t *dst = res->data, *s0, *s1 = e->data, *s2 = f->data;
  element_t *a0, *a1, *b0, *b1;
  element_t p0, p1, p2, p3;

  a0 = s1;
  a1 = &s1[3];
  b0 = s2;
  b1 = &s2[3];

  element_init(p0, res->field);
  element_init(p1, res->field);
  element_init(p2, res->field);
  element_init(p3, res->field);

  s0 = p0->data;
  s1 = p1->data;
  s2 = p2->data;
  element_add(s0[0], a0[0], a1[0]);
  element_add(s0[1], a0[1], a1[1]);
  element_add(s0[2], a0[2], a1[2]);

  element_add(s1[0], b0[0], b1[0]);
  element_add(s1[1], b0[1], b1[1]);
  element_add(s1[2], b0[2], b1[2]);

  kar_poly_2(s2, s2[3], s2[4], s0, s1, p3->data);
  kar_poly_2(s0, s0[3], s0[4], a0, b0, p3->data);
  kar_poly_2(s1, s1[3], s1[4], a1, b1, p3->data);

  element_set(dst[0], s0[0]);
  element_set(dst[1], s0[1]);
  element_set(dst[2], s0[2]);

  element_sub(dst[3], s0[3], s0[0]);
  element_sub(dst[3], dst[3], s1[0]);
  element_add(dst[3], dst[3], s2[0]);

  element_sub(dst[4], s0[4], s0[1]);
  element_sub(dst[4], dst[4], s1[1]);
  element_add(dst[4], dst[4], s2[1]);

  element_sub(dst[5], s2[2], s0[2]);
  element_sub(dst[5], dst[5], s1[2]);

  // Start reusing part of s0 as scratch space(!)
  element_sub(s0[0], s2[3], s0[3]);
  element_sub(s0[0], s0[0], s1[3]);
  element_add(s0[0], s0[0], s1[0]);

  element_sub(s0[1], s2[4], s0[4]);
  element_sub(s0[1], s0[1], s1[4]);
  element_add(s0[1], s0[1], s1[1]);

  polymod_const_mul(p3, s0[0], p->xpwr[0]);
  element_add(res, res, p3);
  polymod_const_mul(p3, s0[1], p->xpwr[1]);
  element_add(res, res, p3);
  polymod_const_mul(p3, s1[2], p->xpwr[2]);
  element_add(res, res, p3);
  polymod_const_mul(p3, s1[3], p->xpwr[3]);
  element_add(res, res, p3);
  polymod_const_mul(p3, s1[4], p->xpwr[4]);
  element_add(res, res, p3);

  element_clear(p0);
  element_clear(p1);
  element_clear(p2);
  element_clear(p3);
}

// General polynomial modulo ring multiplication.
static void polymod_mul(element_ptr res, element_ptr e, element_ptr f) {
  mfptr p = res->field->data;
  int n = p->n;
  element_t *dst;
  element_t *s1 = e->data, *s2 = f->data;
  element_t prod, p0, c0;
  int i, j;
  element_t *high;  // Coefficients of x^n, ..., x^{2n-2}.

  high = pbc_malloc(sizeof(element_t) * (n - 1));
  for (i=0; i<n-1; i++) {
    element_init(high[i], p->field);
    element_set0(high[i]);
  }
  element_init(prod, res->field);
  dst = prod->data;
  element_init(p0, res->field);
  element_init(c0, p->field);

  for (i=0; i<n; i++) {
    int ni = n - i;
    for (j=0; j<ni; j++) {
      element_mul(c0, s1[i], s2[j]);
      element_add(dst[i + j], dst[i + j], c0);
    }
    for (;j<n; j++) {
      element_mul(c0, s1[i], s2[j]);
      element_add(high[j - ni], high[j - ni], c0);
    }
  }

  for (i=0; i<n-1; i++) {
    polymod_const_mul(p0, high[i], p->xpwr[i]);
    element_add(prod, prod, p0);
    element_clear(high[i]);
  }
  pbc_free(high);

  element_set(res, prod);
  element_clear(prod);
  element_clear(p0);
  element_clear(c0);
}

static void polymod_square_degree3(element_ptr res, element_ptr e) {
  // TODO: Investigate if squaring is significantly cheaper than
  // multiplication. If so convert to Karatsuba.
  element_t *dst = res->data;
  element_t *src = e->data;
  mfptr p = res->field->data;
  element_t p0;
  element_t c0, c2;
  element_ptr c1, c3;

  element_init(p0, res->field);
  element_init(c0, p->field);
  element_init(c2, p->field);

  c3 = p0->data;
  c1 = c3 + 1;

  element_mul(c3, src[0], src[1]);
  element_mul(c1, src[0], src[2]);
  element_square(dst[0], src[0]);

  element_mul(c2, src[1], src[2]);
  element_square(c0, src[2]);
  element_square(dst[2], src[1]);

  element_add(dst[1], c3, c3);

  element_add(c1, c1, c1);
  element_add(dst[2], dst[2], c1);

  polymod_const_mul(p0, c0, p->xpwr[1]);
  element_add(res, res, p0);

  element_add(c2, c2, c2);
  polymod_const_mul(p0, c2, p->xpwr[0]);
  element_add(res, res, p0);

  element_clear(p0);
  element_clear(c0);
  element_clear(c2);
}

static void polymod_square(element_ptr res, element_ptr e) {
  element_t *dst;
  element_t *src = e->data;
  mfptr p = res->field->data;
  int n = p->n;
  element_t prod, p0, c0;
  int i, j;
  element_t *high; // Coefficients of x^n,...,x^{2n-2}.

  high = pbc_malloc(sizeof(element_t) * (n - 1));
  for (i=0; i<n-1; i++) {
    element_init(high[i], p->field);
    element_set0(high[i]);
  }

  element_init(prod, res->field);
  dst = prod->data;
  element_init(p0, res->field);
  element_init(c0, p->field);

  for (i=0; i<n; i++) {
    int twicei = 2 * i;
    element_square(c0, src[i]);
    if (twicei < n) {
      element_add(dst[twicei], dst[twicei], c0);
    } else {
      element_add(high[twicei - n], high[twicei - n], c0);
    }

    for (j=i+1; j<n-i; j++) {
      element_mul(c0, src[i], src[j]);
      element_add(c0, c0, c0);
      element_add(dst[i + j], dst[i + j], c0);
    }
    for (;j<n; j++) {
      element_mul(c0, src[i], src[j]);
      element_add(c0, c0, c0);
      element_add(high[i + j - n], high[i + j - n], c0);
    }
  }

  for (i=0; i<n-1; i++) {
    polymod_const_mul(p0, high[i], p->xpwr[i]);
    element_add(prod, prod, p0);
    element_clear(high[i]);
  }
  pbc_free(high);

  element_set(res, prod);
  element_clear(prod);
  element_clear(p0);
  element_clear(c0);
}

static int polymod_is0(element_ptr e) {
  mfptr p = e->field->data;
  element_t *coeff = e->data;
  int i, n = p->n;

  for (i=0; i<n; i++) {
    if (!element_is0(coeff[i])) return 0;
  }
  return 1;
}

static int polymod_is1(element_ptr e) {
  mfptr p = e->field->data;
  element_t *coeff = e->data;
  int i, n = p->n;

  if (!element_is1(coeff[0])) return 0;
  for (i=1; i<n; i++) {
    if (!element_is0(coeff[i])) return 0;
  }
  return 1;
}

static void polymod_set0(element_ptr e) {
  mfptr p = e->field->data;
  element_t *coeff = e->data;
  int i, n = p->n;

  for (i=0; i<n; i++) {
    element_set0(coeff[i]);
  }
}

static void polymod_set1(element_ptr e) {
  mfptr p = e->field->data;
  element_t *coeff = e->data;
  int i, n = p->n;

  element_set1(coeff[0]);
  for (i=1; i<n; i++) {
    element_set0(coeff[i]);
  }
}

static int polymod_sgn(element_ptr e) {
  mfptr p = e->field->data;
  element_t *coeff = e->data;
  int res = 0;
  int i, n = p->n;
  for (i=0; i<n; i++) {
    res = element_sgn(coeff[i]);
    if (res) break;
  }
  return res;
}

static size_t polymod_out_str(FILE *stream, int base, element_ptr e) {
  size_t result = 2, status;
  mfptr p = e->field->data;
  element_t *coeff = e->data;
  int i, n = p->n;

  if (EOF == fputc('[', stream)) return 0;
  for (i=0; i<n; i++) {
    if (i) {
      if (EOF == fputs(", ", stream)) return 0;
      result += 2;
    }
    status = element_out_str(stream, base, coeff[i]);
    if (!status) return 0;
    result += status;
  }
  if (EOF == fputc(']', stream)) return 0;
  return result;
}

static int polymod_snprint(char *s, size_t size, element_ptr e) {
  mfptr p = e->field->data;
  element_t *coeff = e->data;
  int i, n = p->n;
  size_t result = 0, left;
  int status;

  #define clip_sub(void) {                     \
    result += status;                          \
    left = result >= size ? 0 : size - result; \
  }

  status = snprintf(s, size, "[");
  if (status < 0) return status;
  clip_sub();

  for (i=0; i<n; i++) {
    if (i) {
      status = snprintf(s + result, left, ", ");
      if (status < 0) return status;
      clip_sub();
    }
    status = element_snprint(s + result, left, coeff[i]);
    if (status < 0) return status;
    clip_sub();
  }
  status = snprintf(s + result, left, "]");
  if (status < 0) return status;
  return result + status;
  #undef clip_sub
}

static void polymod_set_multiz(element_ptr e, multiz m) {
  mfptr p = e->field->data;
  element_t *coeff = e->data;
  int i, n = p->n;
  if (multiz_is_z(m)) {
    element_set_multiz(coeff[0], m);
    for (i = 1; i < n; i++) element_set0(coeff[i]);
    return;
  }
  int max = multiz_count(m);
  for (i = 0; i < n; i++) {
    if (i >= max) element_set0(coeff[i]);
    else element_set_multiz(coeff[i], multiz_at(m, i));
  }
}

static int polymod_set_str(element_ptr e, const char *s, int base) {
  mfptr p = e->field->data;
  element_t *coeff = e->data;
  int i, n = p->n;
  const char *cp = s;
  element_set0(e);
  while (*cp && isspace(*cp)) cp++;
  if (*cp++ != '[') return 0;
  for (i=0; i<n; i++) {
    cp += element_set_str(coeff[i], cp, base);
    while (*cp && isspace(*cp)) cp++;
    if (i<n-1 && *cp++ != ',') return 0;
  }
  if (*cp++ != ']') return 0;
  return cp - s;
}

static int polymod_coeff_count(element_ptr e) {
  UNUSED_VAR(e);
  mfptr p = e->field->data;
  return p->n;
}

static element_ptr polymod_coeff(element_ptr e, int i) {
  element_t *coeff = e->data;
  return coeff[i];
}

static void polymod_to_mpz(mpz_t z, element_ptr e) {
  element_to_mpz(z, polymod_coeff(e, 0));
}

// Compute x^n,...,x^{2n-2} mod poly.
static void compute_x_powers(field_ptr field, element_ptr poly) {
  mfptr p = field->data;
  element_t p0;
  element_ptr pwrn;
  element_t *coeff, *coeff1;
  int i, j;
  int n = p->n;
  element_t *xpwr;

  xpwr = p->xpwr;

  element_init(p0, field);
  for (i=0; i<n; i++) {
    element_init(xpwr[i], field);
  }
  pwrn = xpwr[0];
  poly_to_polymod_truncate(pwrn, poly);
  element_neg(pwrn, pwrn);

  for (i=1; i<n; i++) {
    coeff = xpwr[i-1]->data;
    coeff1 = xpwr[i]->data;

    element_set0(coeff1[0]);
    for (j=1; j<n; j++) {
      element_set(coeff1[j], coeff[j - 1]);
    }
    polymod_const_mul(p0, coeff[n - 1], pwrn);
    element_add(xpwr[i], xpwr[i], p0);
  }
  element_clear(p0);
}

static void polymod_out_info(FILE *str, field_ptr f) {
  mfptr p = f->data;
  element_fprintf(str, "Extension, poly = %B, base field = ", p->poly);
  field_out_info(str, p->field);
}

// Sets d = gcd(f, g).
static void poly_gcd(element_ptr d, element_ptr f, element_ptr g) {
  element_t a, b, q, r;
  element_init(a, d->field);
  element_init(b, d->field);
  element_init(q, d->field);
  element_init(r, d->field);

  element_set(a, f);
  element_set(b, g);
  for(;;) {
    //TODO: don't care about q
    poly_div(q, r, a, b);
    if (element_is0(r)) break;
    element_set(a, b);
    element_set(b, r);
  }
  element_set(d, b);
  element_clear(a);
  element_clear(b);
  element_clear(q);
  element_clear(r);
}

// Sets f = c g where c is the inverse of the leading coefficient of g.
static void poly_make_monic(element_t f, element_t g) {
  int n = poly_coeff_count(g);
  int i;
  element_ptr e0;
  poly_alloc(f, n);
  if (!n) return;

  e0 = poly_coeff(f, n - 1);
  element_invert(e0, poly_coeff(g, n - 1));
  for (i=0; i<n-1; i++) {
    element_mul(poly_coeff(f, i), poly_coeff(g, i), e0);
  }
  element_set1(e0);
}

// The above should be static.

void field_init_poly(field_ptr f, field_ptr base_field) {
  field_init(f);
  pfptr p = f->data = pbc_malloc(sizeof(*p));
  p->field = base_field;
  p->mapbase = element_field_to_poly;
  f->field_clear = field_clear_poly;
  f->init = poly_init;
  f->clear = poly_clear;
  f->set_si = poly_set_si;
  f->set_multiz = poly_set_multiz;
  f->set_mpz = poly_set_mpz;
  f->to_mpz = poly_to_mpz;
  f->out_str = poly_out_str;
  f->snprint = poly_snprint;
  f->set = poly_set;
  f->sign = poly_sgn;
  f->add = poly_add;
  f->doub = poly_double;
  f->is0 = poly_is0;
  f->is1 = poly_is1;
  f->set0 = poly_set0;
  f->set1 = poly_set1;
  f->sub = poly_sub;
  f->neg = poly_neg;
  f->mul = poly_mul;
  f->mul_mpz = poly_mul_mpz;
  f->mul_si = poly_mul_si;
  f->cmp = poly_cmp;
  f->out_info = poly_out_info;
  f->item_count = poly_coeff_count;
  f->item = poly_coeff;

  f->to_bytes = poly_to_bytes;
  f->from_bytes = poly_from_bytes;
  f->fixed_length_in_bytes = -1;
  f->length_in_bytes = poly_length_in_bytes;
}

void poly_set_coeff(element_ptr e, element_ptr a, int n) {
  peptr p = e->data;
  if (p->coeff->count < n + 1) {
    poly_alloc(e, n + 1);
  }
  element_ptr e0 = p->coeff->item[n];
  element_set(e0, a);
  if (p->coeff->count == n + 1 && element_is0(a)) poly_remove_leading_zeroes(e);
}

void poly_set_coeff0(element_ptr e, int n) {
  peptr p = e->data;
  if (n < p->coeff->count) {
    element_set0(p->coeff->item[n]);
    if (n == p->coeff->count - 1) poly_remove_leading_zeroes(e);
  }
}

void poly_set_coeff1(element_ptr e, int n) {
  peptr p = e->data;
  if (p->coeff->count < n + 1) {
    poly_alloc(e, n + 1);
  }
  element_set1(p->coeff->item[n]);
}

void poly_setx(element_ptr f) {
  poly_alloc(f, 2);
  element_set1(poly_coeff(f, 1));
  element_set0(poly_coeff(f, 0));
}

void poly_const_mul(element_ptr res, element_ptr a, element_ptr poly) {
  int i, n = poly_coeff_count(poly);
  poly_alloc(res, n);
  for (i=0; i<n; i++) {
    element_mul(poly_coeff(res, i), a, poly_coeff(poly, i));
  }
  poly_remove_leading_zeroes(res);
}

void poly_random_monic(element_ptr f, int deg) {
  int i;
  poly_alloc(f, deg + 1);
  for (i=0; i<deg; i++) {
    element_random(poly_coeff(f, i));
  }
  element_set1(poly_coeff(f, i));
}

int polymod_field_degree(field_t f) {
  mfptr p = f->data;
  return p->n;
}

void field_init_polymod(field_ptr f, element_ptr poly) {
  pfptr pdp = poly->field->data;
  field_init(f);
  mfptr p = f->data = pbc_malloc(sizeof(*p));
  p->field = pdp->field;
  p->mapbase = element_field_to_poly;
  element_init(p->poly, poly->field);
  element_set(p->poly, poly);
  int n = p->n = poly_degree(p->poly);
  f->field_clear = field_clear_polymod;
  f->init = polymod_init;
  f->clear = polymod_clear;
  f->set_si = polymod_set_si;
  f->set_mpz = polymod_set_mpz;
  f->out_str = polymod_out_str;
  f->snprint = polymod_snprint;
  f->set_multiz = polymod_set_multiz;
  f->set_str = polymod_set_str;
  f->set = polymod_set;
  f->sign = polymod_sgn;
  f->add = polymod_add;
  f->doub = polymod_double;
  f->sub = polymod_sub;
  f->neg = polymod_neg;
  f->is0 = polymod_is0;
  f->is1 = polymod_is1;
  f->set0 = polymod_set0;
  f->set1 = polymod_set1;
  f->cmp = polymod_cmp;
  f->to_mpz = polymod_to_mpz;
  f->item_count = polymod_coeff_count;
  f->item = polymod_coeff;
  switch(n) {
    case 3:
      f->mul = polymod_mul_degree3;
      f->square = polymod_square_degree3;
      break;
    case 6:
      f->mul = polymod_mul_degree6;
      f->square = polymod_square;
      break;
    default:
      f->mul = polymod_mul;
      f->square = polymod_square;
      break;
  }

  f->mul_mpz = polymod_mul_mpz;
  f->mul_si = polymod_mul_si;
  f->random = polymod_random;
  f->from_hash = polymod_from_hash;
  f->invert = polymod_invert;
  f->is_sqr = polymod_is_sqr;
  f->sqrt = polymod_sqrt;
  f->to_bytes = polymod_to_bytes;
  f->from_bytes = polymod_from_bytes;
  f->out_info = polymod_out_info;

  if (pdp->field->fixed_length_in_bytes < 0) {
    f->fixed_length_in_bytes = -1;
    f->length_in_bytes = polymod_length_in_bytes;
  } else {
    f->fixed_length_in_bytes = pdp->field->fixed_length_in_bytes * poly_degree(poly);
  }
  mpz_pow_ui(f->order, p->field->order, n);

  p->xpwr = pbc_malloc(sizeof(element_t) * n);
  compute_x_powers(f, poly);
}

field_ptr poly_base_field(element_t f) {
  return ((pfptr) f->field->data)->field;
}

void polymod_const_mul(element_ptr res, element_ptr a, element_ptr e) {
  // a lies in R, e in R[x].
  element_t *coeff = e->data, *dst = res->data;
  int i, n = polymod_field_degree(e->field);

  for (i=0; i<n; i++) {
    element_mul(dst[i], coeff[i], a);
  }
}

struct checkgcd_scope_var {
  mpz_ptr     z, deg;
  field_ptr   basef;
  element_ptr xpow, x, f, g;
};

// Returns 0 if gcd(x^q^{n/d} - x, f) = 1, 1 otherwise.
static int checkgcd(mpz_ptr fac, unsigned int mul, struct checkgcd_scope_var *v) {
  UNUSED_VAR(mul);
  mpz_divexact(v->z, v->deg, fac);
  mpz_pow_ui(v->z, v->basef->order, mpz_get_ui(v->z));
  element_pow_mpz(v->xpow, v->x, v->z);
  element_sub(v->xpow, v->xpow, v->x);
  if (element_is0(v->xpow)) return 1;
  polymod_to_poly(v->g, v->xpow);
  poly_gcd(v->g, v->f, v->g);
  return poly_degree(v->g) != 0;
}

// Returns 1 if polynomial is irreducible, 0 otherwise.
// A polynomial f(x) is irreducible in F_q[x] if and only if:
//  (1) f(x) | x^{q^n} - x, and
//  (2) gcd(f(x), x^{q^{n/d}} - x) = 1 for all primes d | n.
// (Recall GF(p) is the splitting field for x^p - x.)
int poly_is_irred(element_ptr f) {
  int res = 0;
  element_t xpow, x, g;
  field_ptr basef = poly_base_field(f);
  field_t rxmod;

  // 0, units are not irreducibles.
  // Assume coefficients are from a field.
  if (poly_degree(f) <= 0) return 0;
  // Degree 1 polynomials are always irreducible.
  if (poly_degree(f) == 1) return 1;

  field_init_polymod(rxmod, f);
  element_init(xpow, rxmod);
  element_init(x, rxmod);
  element_init(g, f->field);
  element_set1(polymod_coeff(x, 1));

  // The degree fits in an unsigned int but I'm lazy and want to use my
  // mpz trial division code.
  mpz_t deg, z;
  mpz_init(deg);
  mpz_init(z);
  mpz_set_ui(deg, poly_degree(f));

  struct checkgcd_scope_var v = {.z = z, .deg = deg, .basef = basef,
                                 .xpow = xpow, .x = x, .f = f, .g = g};
  if (!pbc_trial_divide((int(*)(mpz_t,unsigned,void*))checkgcd, &v, deg, NULL)) {
    // By now condition (2) has been satisfied. Check (1).
    mpz_pow_ui(z, basef->order, poly_degree(f));
    element_pow_mpz(xpow, x, z);
    element_sub(xpow, xpow, x);
    if (element_is0(xpow)) res = 1;
  }

  mpz_clear(deg);
  mpz_clear(z);
  element_clear(g);
  element_clear(xpow);
  element_clear(x);
  field_clear(rxmod);
  return res;
}

void element_field_to_poly(element_ptr f, element_ptr g) {
  poly_alloc(f, 1);
  element_set(poly_coeff(f, 0), g);
  poly_remove_leading_zeroes(f);
}

void element_field_to_polymod(element_ptr f, element_ptr g) {
  mfptr p = f->field->data;
  element_t *coeff = f->data;
  int i, n = p->n;
  element_set(coeff[0], g);
  for (i=1; i<n; i++) {
    element_set0(coeff[i]);
  }
}

// Returns 0 when a root exists and sets root to one of the roots.
int poly_findroot(element_ptr root, element_ptr poly) {
  // Compute gcd(x^q - x, poly).
  field_t fpxmod;
  element_t p, x, r, fac, g;
  mpz_t q;

  mpz_init(q);
  mpz_set(q, poly_base_field(poly)->order);

  field_init_polymod(fpxmod, poly);
  element_init(p, fpxmod);
  element_init(x, fpxmod);
  element_init(g, poly->field);
  element_set1(((element_t *) x->data)[1]);
pbc_info("findroot: degree %d...", poly_degree(poly));
  element_pow_mpz(p, x, q);
  element_sub(p, p, x);

  polymod_to_poly(g, p);
  element_clear(p);
  poly_gcd(g, g, poly);
  poly_make_monic(g, g);
  element_clear(x);
  field_clear(fpxmod);

  if (!poly_degree(g)) {
    printf("no roots!\n");
    mpz_clear(q);
    element_clear(g);
    return -1;
  }

  // Cantor-Zassenhaus algorithm.
  element_init(fac, g->field);
  element_init(x, g->field);
  element_set_si(x, 1);
  mpz_sub_ui(q, q, 1);
  mpz_divexact_ui(q, q, 2);
  element_init(r, g->field);
  for (;;) {
    if (poly_degree(g) == 1) break;  // Found a root!
step_random:
    poly_random_monic(r, 1);
    // TODO: evaluate at g instead of bothering with gcd
    poly_gcd(fac, r, g);

    if (poly_degree(fac) > 0) {
      poly_make_monic(g, fac);
    } else {
      field_init_polymod(fpxmod, g);
      int n;
      element_init(p, fpxmod);

      poly_to_polymod_truncate(p, r);
pbc_info("findroot: degree %d...", poly_degree(g));
      element_pow_mpz(p, p, q);

      polymod_to_poly(r, p);
      element_clear(p);
      field_clear(fpxmod);

      element_add(r, r, x);
      poly_gcd(fac, r, g);
      n = poly_degree(fac);
      if (n > 0 && n < poly_degree(g)) {
        poly_make_monic(g, fac);
      } else {
        goto step_random;
      }
    }
  }
pbc_info("findroot: found root");
  element_neg(root, poly_coeff(g, 0));
  element_clear(r);
  mpz_clear(q);
  element_clear(x);
  element_clear(g);
  element_clear(fac);
  return 0;
}
