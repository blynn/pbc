// Test polynomials.
#include "pbc.h"
#include "pbc_fp.h"
#include "pbc_poly.h"
#include "pbc_test.h"
#include "misc/darray.h"

static void elfree(void *data) {
  element_clear(data);
  pbc_free(data);
}

static void inner(void *data2, element_ptr f, field_t fx, darray_t prodlist) {
  element_ptr g = data2;
  if (!poly_degree(f) || !poly_degree(g)) return;
  if (poly_degree(f) + poly_degree(g) > 3) return;
  element_ptr h = pbc_malloc(sizeof(*h));
  element_init(h, fx);
  element_mul(h, f, g);
  darray_append(prodlist, h);
  EXPECT(!poly_is_irred(h));
}

static void outer(void *data, darray_t list, field_t fx, darray_t prodlist) {
  element_ptr f = data;
  darray_forall4(list, (void(*)(void*,void*,void*,void*))inner, f, fx, prodlist);
}

int isf(void *data, element_ptr f) {
  element_ptr f1 = data;
  return !element_cmp(f, f1);
}

int main(void) {
  field_t fp, fx;
  mpz_t prime;
  darray_t list;
  int p = 7;

  // Exercise poly_is_irred() with a sieve of Erastosthenes for polynomials.
  darray_init(list);
  mpz_init(prime);
  mpz_set_ui(prime, p);
  field_init_fp(fp, prime);
  field_init_poly(fx, fp);
  element_t e;
  element_init(e, fp);
  // Enumerate polynomials in F_p[x] up to degree 2.
  int a[3], d;
  a[0] = a[1] = a[2] = 0;
  for(;;) {
    element_ptr f = pbc_malloc(sizeof(*f));
    element_init(f, fx);
    int j;
    for(j = 0; j < 3; j++) {
      element_set_si(e, a[j]);
      poly_set_coeff(f, e, j);
    }

    // Test poly_degree().
    for(j = 2; j >= 0 && !a[j]; j--);
    EXPECT(poly_degree(f) == j);

    // Add monic polynomials to the list.
    if (j >= 0 && a[j] == 1) darray_append(list, f);
    else {
      element_clear(f);
      pbc_free(f);
    }

    // Next!
    d = 0;
    for(;;) {
      a[d]++;
      if (a[d] >= p) {
        a[d] = 0;
        d++;
        if (d > 2) goto break2;
      } else break;
    }
  }
break2: ;

  // Find all composite monic polynomials of degree 3 or less.
  darray_t prodlist;
  darray_init(prodlist);

  darray_forall4(list, (void(*)(void*,void*,void*,void*))outer, list, fx, prodlist);

  // Enumerate all monic polynomials in F_p[x] up to degree 3.
  a[0] = a[1] = a[2] = 0;
  for(;;) {
    element_t f;
    element_init(f, fx);
    int j;
    for(j = 0; j < 3; j++) {
      element_set_si(e, a[j]);
      poly_set_coeff(f, e, j);
    }
    for(j = 2; j >= 0 && !a[j]; j--);
    element_set1(e);
    poly_set_coeff(f, e, j + 1);

    // Check f is a unit or appears on the list of composites if and only if
    // poly_is_irred() returns 0.
    if (poly_is_irred(f)) {
      EXPECT(!darray_at_test(prodlist, (int(*)(void*,void*))isf, f));
    } else if (poly_degree(f)) {
      EXPECT(darray_at_test(prodlist, (int(*)(void*,void*))isf, f));
    }
    element_clear(f);

    // Next!
    d = 0;
    for(;;) {
      a[d]++;
      if (a[d] >= p) {
        a[d] = 0;
        d++;
        if (d > 2) goto break3;
      } else break;
    }
  }
break3: ;

  darray_forall(list, elfree);
  darray_forall(prodlist, elfree);
  darray_clear(prodlist);
  darray_clear(list);
  mpz_clear(prime);
  field_clear(fx);
  field_clear(fp);
  element_clear(e);

  return pbc_err_count;
}
