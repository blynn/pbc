// Routines for finding:
// * MNT curves with embedding degree 6
// * Freeman curves (which have embedding degree 10)

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h> // for intptr_t
#include <gmp.h>
#include "pbc_mnt.h"
#include "pbc_memory.h"
#include "pbc_utils.h"
#include "misc/darray.h"

struct pell_solution_s {
  int count;
  mpz_t minx; //minimal solution of x^2 - Dy^2 = 1
  mpz_t miny;
  mpz_t *x;
  mpz_t *y;
};
typedef struct pell_solution_s pell_solution_t[1];
typedef struct pell_solution_s *pell_solution_ptr;

static void freempz(void *data) {
  mpz_clear(data);
  pbc_free(data);
}

// Solves x^2 - Dy^2 = N where D not a square.
// For square D, we have (x+Dy)(x-Dy) = N so we look at the factors of N.
static void general_pell(pell_solution_t ps, mpz_t D, int N) {
  // TODO: Use brute force for small D.
  int i, sgnN = N > 0 ? 1 : -1;
  intptr_t f, n;

  // Find square factors of N.
  darray_t listf;
  darray_init(listf);

  f = 1;
  for (;;) {
    n = f * f;
    if (n > abs(N)) break;
    if (!(abs(N) % n)) {
      darray_append(listf, int_to_voidp(f));
    }
    f++;
  }

  //a0, twice_a0 don't change once initialized
  //a1 is a_i every iteration
  //P0, P1 become P_{i-1}, P_i every iteration
  //similarly for Q0, Q1
  mpz_t a0, twice_a0, a1;
  mpz_t P0, P1;
  mpz_t Q0, Q1;
  //variables to compute the convergents
  mpz_t p0, p1, pnext;
  mpz_t q0, q1, qnext;

  int d;

  darray_t listp, listq;
  mpz_ptr zptr;

  mpz_init(a0);
  mpz_init(twice_a0);
  mpz_init(a1);
  mpz_init(P0); mpz_init(P1);
  mpz_init(Q0); mpz_init(Q1);
  mpz_init(p0); mpz_init(p1); mpz_init(pnext);
  mpz_init(q0); mpz_init(q1); mpz_init(qnext);

  darray_init(listp);
  darray_init(listq);

  mpz_sqrt(a0, D);
  mpz_set_ui(P0, 0);
  mpz_set_ui(Q0, 1);

  mpz_set(P1, a0);
  mpz_mul(Q1, a0, a0);
  mpz_sub(Q1, D, Q1);
  mpz_add(a1, a0, P1);
  mpz_tdiv_q(a1, a1, Q1);

  mpz_add(twice_a0, a0, a0);

  mpz_set(p0, a0);
  mpz_set_ui(q0, 1);
  mpz_mul(p1, a0, a1);
  mpz_add_ui(p1, p1, 1);
  mpz_set(q1, a1);

  d = -1;
  for(;;) {
    if (d == sgnN) {
      for (i=0; i<listf->count; i++) {
        f = (intptr_t) listf->item[i];
        if (!mpz_cmp_ui(Q1, abs(N) / (f * f))) {
//element_printf("found %Zd, %Zd, %d\n", p0, q0, f);
          zptr = (mpz_ptr) pbc_malloc(sizeof(mpz_t));
          mpz_init(zptr);
          mpz_set(zptr, p0);
          mpz_mul_ui(zptr, p0, f);
          darray_append(listp, zptr);
          zptr = (mpz_ptr) pbc_malloc(sizeof(mpz_t));
          mpz_init(zptr);
          mpz_set(zptr, q0);
          mpz_mul_ui(zptr, q0, f);
          darray_append(listq, zptr);
        }
      }
    }

    if (!mpz_cmp(twice_a0, a1) && d == 1) break;
    //compute more of the continued fraction expansion
    mpz_set(P0, P1);
    mpz_mul(P1, a1, Q1);
    mpz_sub(P1, P1, P0);
    mpz_set(Q0, Q1);
    mpz_mul(Q1, P1, P1);
    mpz_sub(Q1, D, Q1);
    mpz_divexact(Q1, Q1, Q0);
    mpz_add(a1, a0, P1);
    mpz_tdiv_q(a1, a1, Q1);

    //compute next convergent
    mpz_mul(pnext, a1, p1);
    mpz_add(pnext, pnext, p0);
    mpz_set(p0, p1);
    mpz_set(p1, pnext);

    mpz_mul(qnext, a1, q1);
    mpz_add(qnext, qnext, q0);
    mpz_set(q0, q1);
    mpz_set(q1, qnext);
    d = -d;
  }
  darray_clear(listf);

  mpz_init(ps->minx);
  mpz_init(ps->miny);
  mpz_set(ps->minx, p0);
  mpz_set(ps->miny, q0);
  n = listp->count;
  ps->count = n;
  if (n) {
    ps->x = (mpz_t *) pbc_malloc(sizeof(mpz_t) * n);
    ps->y = (mpz_t *) pbc_malloc(sizeof(mpz_t) * n);
    for (i = 0; i < n; i++) {
      mpz_init(ps->x[i]);
      mpz_init(ps->y[i]);
      mpz_set(ps->x[i], (mpz_ptr) listp->item[i]);
      mpz_set(ps->y[i], (mpz_ptr) listq->item[i]);
    }
  }

  mpz_clear(a0);
  mpz_clear(twice_a0);
  mpz_clear(a1);
  mpz_clear(P0); mpz_clear(P1);
  mpz_clear(Q0); mpz_clear(Q1);
  mpz_clear(p0); mpz_clear(p1); mpz_clear(pnext);
  mpz_clear(q0); mpz_clear(q1); mpz_clear(qnext);

  darray_forall(listp, freempz);
  darray_forall(listq, freempz);
  darray_clear(listp);
  darray_clear(listq);
}

static void pell_solution_clear(pell_solution_t ps) {
  int i, n = ps->count;

  if (n) {
    for (i=0; i<n; i++) {
      mpz_clear(ps->x[i]);
      mpz_clear(ps->y[i]);
    }
    pbc_free(ps->x);
    pbc_free(ps->y);
  }
  mpz_clear(ps->minx);
  mpz_clear(ps->miny);
}

void pbc_cm_init(pbc_cm_t cm) {
  mpz_init(cm->q);
  mpz_init(cm->r);
  mpz_init(cm->h);
  mpz_init(cm->n);
}

void pbc_cm_clear(pbc_cm_t cm) {
  mpz_clear(cm->q);
  mpz_clear(cm->r);
  mpz_clear(cm->h);
  mpz_clear(cm->n);
}

static int mnt_step2(int (*callback)(pbc_cm_t, void *), void *data,
    unsigned int D, mpz_t U) {
  int d;
  mpz_t n, l, q;
  mpz_t p;
  mpz_t r, cofac;

  mpz_init(l);
  mpz_mod_ui(l, U, 6);
  if (!mpz_cmp_ui(l, 1)) {
    mpz_sub_ui(l, U, 1);
    d = 1;
  } else if (!mpz_cmp_ui(l, 5)) {
    mpz_add_ui(l, U, 1);
    d = -1;
  } else {
    mpz_clear(l);
    return 0;
  }

  mpz_divexact_ui(l, l, 3);
  mpz_init(q);

  mpz_mul(q, l, l);
  mpz_add_ui(q, q, 1);
  if (!mpz_probab_prime_p(q, 10)) {
    mpz_clear(q);
    mpz_clear(l);
    return 0;
  }

  mpz_init(n);
  if (d < 0) {
    mpz_sub(n, q, l);
  } else {
    mpz_add(n, q, l);
  }

  mpz_init(p);
  mpz_init(r);
  mpz_init(cofac);
  {
  mpz_set_ui(cofac, 1);
  mpz_set(r, n);
  mpz_set_ui(p, 2);
  if (!mpz_probab_prime_p(r, 10)) for(;;) {
    if (mpz_divisible_p(r, p)) do {
      mpz_mul(cofac, cofac, p);
      mpz_divexact(r, r, p);
    } while (mpz_divisible_p(r, p));
    if (mpz_probab_prime_p(r, 10)) break;
    //TODO: use a table of primes instead?
    mpz_nextprime(p, p);
    if (mpz_sizeinbase(p, 2) > 16) {
      //printf("has 16+ bit factor\n");
      mpz_clear(r);
      mpz_clear(p);
      mpz_clear(cofac);
      mpz_clear(q);
      mpz_clear(l);
      mpz_clear(n);
      return 0;
    }
  }
  }

  pbc_cm_t cm;
  pbc_cm_init(cm);
  cm->k = 6;
  cm->D = D;
  mpz_set(cm->q, q);
  mpz_set(cm->r, r);
  mpz_set(cm->h, cofac);
  mpz_set(cm->n, n);
  int res = callback(cm, data);
  pbc_cm_clear(cm);

  mpz_clear(cofac);
  mpz_clear(r);
  mpz_clear(p);
  mpz_clear(q);
  mpz_clear(l);
  mpz_clear(n);
  return res;
}

int pbc_cm_search_d(int (*callback)(pbc_cm_t, void *), void *data,
    unsigned int D, unsigned int bitlimit) {
  mpz_t D3;
  mpz_t t0, t1, t2;

  mpz_init(D3);
  mpz_set_ui(D3, D * 3);

  if (mpz_perfect_square_p(D3)) {
    // The only squares that differ by 8 are 1 and 9,
    // which we get if U=V=1, D=3, but then l is not an integer.
    mpz_clear(D3);
    return 0;
  }

  mpz_init(t0);
  mpz_init(t1);
  mpz_init(t2);

  pell_solution_t ps;
  general_pell(ps, D3, -8);

  int i, n;
  int res = 0;
  n = ps->count;
  if (n) for (;;) {
    for (i=0; i<n; i++) {
      //element_printf("%Zd, %Zd\n", ps->x[i], ps->y[i]);
      res = mnt_step2(callback, data, D, ps->x[i]);
      if (res) goto toobig;
      //compute next solution as follows
      //if p, q is current solution
      //compute new solution p', q' via
      //(p + q sqrt{3D})(t + u sqrt{3D}) = p' + q' sqrt(3D)
      //where t, u is min. solution to Pell equation
      mpz_mul(t0, ps->minx, ps->x[i]);
      mpz_mul(t1, ps->miny, ps->y[i]);
      mpz_mul(t1, t1, D3);
      mpz_add(t0, t0, t1);
      if (2 * mpz_sizeinbase(t0, 2) > bitlimit + 10) goto toobig;
      mpz_mul(t2, ps->minx, ps->y[i]);
      mpz_mul(t1, ps->miny, ps->x[i]);
      mpz_add(t2, t2, t1);
      mpz_set(ps->x[i], t0);
      mpz_set(ps->y[i], t2);
    }
  }
toobig:

  pell_solution_clear(ps);
  mpz_clear(t0);
  mpz_clear(t1);
  mpz_clear(t2);
  mpz_clear(D3);
  return res;
}

static int freeman_step2(int (*callback)(pbc_cm_t, void *), void *data,
    unsigned int D, mpz_t U) {
  mpz_t n, x, q;
  mpz_t p;
  mpz_t r, cofac;
  pbc_cm_t cm;

  mpz_init(x);
  mpz_mod_ui(x, U, 15);
  if (!mpz_cmp_ui(x, 5)) {
    mpz_sub_ui(x, U, 5);
  } else if (!mpz_cmp_ui(x, 10)) {
    mpz_add_ui(x, U, 5);
  } else {
    pbc_die("should never reach here");
    mpz_clear(x);
    return 0;
  }

  mpz_divexact_ui(x, x, 15);
  mpz_init(q);
  mpz_init(r);

  //q = 25x^4 + 25x^3 + 25x^2 + 10x + 3
  mpz_mul(r, x, x);
  mpz_add(q, x, x);
  mpz_mul_ui(r, r, 5);
  mpz_add(q, q, r);
  mpz_mul(r, r, x);
  mpz_add(q, q, r);
  mpz_mul(r, r, x);
  mpz_add(q, q, r);
  mpz_mul_ui(q, q, 5);
  mpz_add_ui(q, q, 3);

  if (!mpz_probab_prime_p(q, 10)) {
    mpz_clear(q);
    mpz_clear(r);
    mpz_clear(x);
    return 0;
  }

  //t = 10x^2 + 5x + 3
  //n = q - t + 1
  mpz_init(n);

  mpz_mul_ui(n, x, 5);
  mpz_mul(r, n, x);
  mpz_add(r, r, r);
  mpz_add(n, n, r);
  mpz_sub(n, q, n);
  mpz_sub_ui(n, n, 2);

  mpz_init(p);
  mpz_init(cofac);
  {
  mpz_set_ui(cofac, 1);
  mpz_set(r, n);
  mpz_set_ui(p, 2);
  if (!mpz_probab_prime_p(r, 10)) for(;;) {
    if (mpz_divisible_p(r, p)) do {
      mpz_mul(cofac, cofac, p);
      mpz_divexact(r, r, p);
    } while (mpz_divisible_p(r, p));
    if (mpz_probab_prime_p(r, 10)) break;
    //TODO: use a table of primes instead?
    mpz_nextprime(p, p);
    if (mpz_sizeinbase(p, 2) > 16) {
      //printf("has 16+ bit factor\n");
      mpz_clear(r);
      mpz_clear(p);
      mpz_clear(cofac);
      mpz_clear(q);
      mpz_clear(x);
      mpz_clear(n);
      return 0;
    }
  }
  }

  pbc_cm_init(cm);
  cm->k = 10;
  cm->D = D;
  mpz_set(cm->q, q);
  mpz_set(cm->r, r);
  mpz_set(cm->h, cofac);
  mpz_set(cm->n, n);
  int res = callback(cm, data);
  pbc_cm_clear(cm);

  mpz_clear(cofac);
  mpz_clear(r);
  mpz_clear(p);
  mpz_clear(q);
  mpz_clear(x);
  mpz_clear(n);
  return res;
}

int pbc_cm_search_g(int (*callback)(pbc_cm_t, void *), void *data,
    unsigned int D, unsigned int bitlimit) {
  int res = 0;
  mpz_t D15;
  mpz_t t0, t1, t2;

  mpz_init(D15);
  mpz_set_ui(D15, D);
  mpz_mul_ui(D15, D15, 15);
  if (mpz_perfect_square_p(D15)) {
    mpz_clear(D15);
    return 0;
  }

  mpz_init(t0);
  mpz_init(t1);
  mpz_init(t2);

  pell_solution_t ps;
  general_pell(ps, D15, -20);

  int i, n;
  n = ps->count;
  if (n) for (;;) {
    for (i=0; i<n; i++) {
      res = freeman_step2(callback, data, D, ps->x[i]);
      if (res) goto toobig;
      // Compute next solution as follows:
      // If p, q is current solution
      // then compute new solution p', q' via
      //   (p + q sqrt{15D})(t + u sqrt{15D}) = p' + q' sqrt(15D)
      // where t, u is min. solution to Pell equation
      mpz_mul(t0, ps->minx, ps->x[i]);
      mpz_mul(t1, ps->miny, ps->y[i]);
      mpz_mul(t1, t1, D15);
      mpz_add(t0, t0, t1);
      if (2 * mpz_sizeinbase(t0, 2) > bitlimit + 10) goto toobig;
      mpz_mul(t2, ps->minx, ps->y[i]);
      mpz_mul(t1, ps->miny, ps->x[i]);
      mpz_add(t2, t2, t1);
      mpz_set(ps->x[i], t0);
      mpz_set(ps->y[i], t2);
    }
  }
toobig:

  pell_solution_clear(ps);
  mpz_clear(t0);
  mpz_clear(t1);
  mpz_clear(t2);
  mpz_clear(D15);
  return res;
}
