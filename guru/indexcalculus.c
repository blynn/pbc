#include <stdio.h>
#include <stdlib.h>
#include <stdint.h> // for intptr_t
#include <string.h>
#include <math.h>
#include <gmp.h>
#include "pbc.h"
#include "pbc_utils.h"

struct cell_s {
  int ind;
  mpz_t data;
};
typedef struct cell_s *cell_ptr;

static cell_ptr newcell(void)
{
  cell_ptr res;
  res = pbc_malloc(sizeof(struct cell_s));
  //res->data = pbc_malloc(sizeof(mpz_t));
  //mpz_init(res->data);
  mpz_init(res->data);
  return res;
}

static void delcell(void *p)
{
  cell_ptr cp = p;
  mpz_clear(cp->data);
  pbc_free(p);
}

static int is_gen(mpz_t x, mpz_t q, darray_ptr fac, darray_ptr mul) {
  int result = 1;
  mpz_t z;
  mpz_t q1;
  int i;
  UNUSED_VAR(mul);

  mpz_init(z);
  mpz_init(q1);

  mpz_sub_ui(q1, q, 1);
  for (i=0; i<fac->count; i++) {
    mpz_divexact(z, q1, fac->item[i]);
    mpz_powm(z, x, z, q);
    if (!mpz_cmp_ui(z, 1)) {
      result = 0;
      break;
    }
  }

  mpz_clear(q1);
  mpz_clear(z);
  return result;
}

// Garner's Algorithm.
// See Algorithm 14.71, Handbook of Cryptography.
static void CRT(mpz_t x, mpz_ptr *v, mpz_ptr *m, int t) {
  mpz_t u;
  mpz_t C[t];
  int i, j;

  mpz_init(u);
  for (i=1; i<t; i++) {
    mpz_init(C[i]);
    mpz_set_ui(C[i], 1);
    for (j=0; j<i; j++) {
      mpz_invert(u, m[j], m[i]);
      mpz_mul(C[i], C[i], u);
      mpz_mod(C[i], C[i], m[i]);
    }
  }
  mpz_set(u, v[0]);
  mpz_set(x, u);
  for (i=1; i<t; i++) {
    mpz_sub(u, v[i], x);
    mpz_mul(u, u, C[i]);
    mpz_mod(u, u, m[i]);
    for (j=0; j<i; j++) {
      mpz_mul(u, u, m[j]);
    }
    mpz_add(x, x, u);
  }

  for (i=1; i<t; i++) mpz_clear(C[i]);
  mpz_clear(u);
}

//TODO: http://www.cecm.sfu.ca/CAG/abstracts/aaron27Jan06.pdf
//TODO: don't need to store last element of list in row[i]
//TODO: linked lists might be better than dynamic arrays (avoids memmove())
//TODO: allow holes in the table
//(if drought lasts too long)
void index_calculus_step1(mpz_t *ind, int r, mpz_t g, mpz_t q,
    darray_ptr fac, darray_ptr mul) {
  int count = 0;
  int i, j;
  mpz_t z, z0, z1;
  int relcount;
  unsigned int *prime = pbc_malloc(sizeof(unsigned int) * r);
  int bundlecount = (r - 10 + 19) / 20;
  mpz_t *bundle = pbc_malloc(sizeof(mpz_t) * bundlecount);
  int faci;
  mpz_t k, km;

  cell_ptr *rel = pbc_malloc(sizeof(cell_ptr) * r);
  cell_ptr *relm = pbc_malloc(sizeof(cell_ptr) * r);
  //''matrix'' is actually a list of matrices
  //(we solve over different moduli and combine using CRT)
  //darray_t **matrix = pbc_malloc(sizeof(darray_t *) * fac->count);
  darray_t *matrix[fac->count];
  int minfound[fac->count];

  for (i=0; i<r; i++) {
    rel[i] = newcell();
    relm[i] = newcell();
  }
  for (i=0; i<fac->count; i++) {
    //similarly ''row'' refers to a list of rows
    darray_t *row = pbc_malloc(sizeof(darray_t) * r);
    for (j=0; j<r; j++) {
      darray_init(row[j]);
    }
    matrix[i] = row;
    minfound[i] = 0;
  }

  mpz_init(k);
  mpz_init(km);
  mpz_init(z);
  mpz_init(z1);
  mpz_init(z0);

  printf("building prime table...\n");
  prime[0] = 2;
  mpz_set_ui(z, 2);
  for (i=1; i<r; i++) {
    mpz_nextprime(z, z);
    prime[i] = mpz_get_ui(z);
  }

  for (i=0; i<bundlecount; i++) {
    mpz_init(bundle[i]);
    mpz_set_ui(bundle[i], 1);
    for (j=0; j<20; j++) {
      int jj = 10 + 20 * i + j;
      if (jj >= r) break;
      mpz_mul_ui(bundle[i], bundle[i], prime[jj]);
    }
    element_printf("bundle %d: %Zd\n", i, bundle[i]);
  }
  printf("searching for r-smooth numbers\n");

  mpz_set_ui(z, 1);
  mpz_init(k);
  int try = 0;
  do {
    mpz_mul(z, z, g);
    mpz_mod(z, z, q);
    mpz_add_ui(k, k, 1);

    /*
    pbc_mpz_random(k, q);
    mpz_powm(z, g, k, q);
    */

    try++;

    mpz_set(z1, z);
    relcount = 0;
    for (i=0; i<10; i++) {
      if (i >= r) break;
      j = 0;
      while (mpz_divisible_ui_p(z1, prime[i])) {
        mpz_divexact_ui(z1, z1, prime[i]);
        j++;
      }
      if (j) {
        rel[relcount]->ind = i;
        mpz_set_ui(rel[relcount]->data, j);
        relcount++;
        if (!mpz_cmp_ui(z1, 1)) goto found;
      }
    }
    for (i=0; i<bundlecount; i++) {
      mpz_gcd(z0, bundle[i], z1);
      if (mpz_cmp_ui(z0, 1)) {
        int ii;
        for (ii = 0; ii < 20; ii++) {
          int jj = 10 + i * 20 + ii;
          if (jj >= r) break;
          j = 0;
          while (mpz_divisible_ui_p(z1, prime[jj])) {
            mpz_divexact_ui(z1, z1, prime[jj]);
            j++;
          }
          if (j) {
            rel[relcount]->ind = jj;
            mpz_set_ui(rel[relcount]->data, j);
            relcount++;
            if (!mpz_cmp_ui(z1, 1)) goto found;
          }
        }
      }
    }
    continue;
found:

/*
    printf("found r-smooth number after %d tries\n", try);

    gmp_printf("g^%Zd = %Zd:", k, z);
    for (i=0; i<relcount; i++) {
      gmp_printf(" %u:%Zd", rel[i]->ind, rel[i]->data);
    }
    printf("\n");
*/
    try = 0;

    for (faci=0; faci<fac->count; faci++) {
      darray_t *row = matrix[faci];
      mpz_ptr order = fac->item[faci];
      int relmcount = 0;
      mpz_mod(km, k, order);

      //gmp_printf("mod %Zd\n", order);
      for (i=0; i<relcount; i++) {
        mpz_mod(z0, rel[i]->data, order);
        if (mpz_sgn(z0)) {
          mpz_set(relm[relmcount]->data, z0);
          relm[relmcount]->ind = rel[i]->ind;
          relmcount++;
        }
      }

      while (relmcount) {
        //start from the sparse end
        int rind = relm[relmcount - 1]->ind;
        darray_ptr rp = row[rind];

        if (rind < minfound[faci]) break;

        mpz_set(z0, relm[relmcount - 1]->data);
        if (!rp->count) {
          mpz_invert(z0, z0, order);
          cell_ptr cnew = newcell();
          cnew->ind = -1;
          mpz_mul(z1, km, z0);
          mpz_mod(cnew->data, z1, order);
          darray_append(rp, cnew);
          for (j=0; j<relmcount; j++) {
            cnew = newcell();
            cnew->ind = relm[j]->ind;
            mpz_mul(z1, relm[j]->data, z0);
            mpz_mod(cnew->data, z1, order);
            darray_append(rp, cnew);
          }
          count++;
printf("%d / %d\n", count, r * fac->count);
/*
for (i=1; i<rp->count; i++) {
  cnew = rp->item[i];
  gmp_printf(" %u:%Zd", cnew->ind, cnew->data);
}
cnew = rp->item[0];
gmp_printf(" %Zd\n", cnew->data);
*/

          if (rind == minfound[faci]) {
            do {
              if (!minfound[faci]) {
              printf("found log p_%d\n", minfound[faci]);
              cnew = rp->item[0];
              gmp_printf("km = %Zd mod %Zd\n", cnew->data, order);
              }
              minfound[faci]++;
              if (minfound[faci] >= r) break;
              rp = row[minfound[faci]];
            } while (rp->count);
          }
          break;

        }

/*
{
//gmp_printf("mod = %Zd\n", order);
printf("before:");
for (i=0; i<relmcount; i++) {
  gmp_printf(" %u:%Zd", relm[i]->ind, relm[i]->data);
}
gmp_printf(" %Zd\n", km);
cell_ptr cp;
printf("sub %d:", rind);
for (i=1; i<rp->count; i++) {
  cp = rp->item[i];
  gmp_printf(" %u:%Zd", cp->ind, cp->data);
}
cp = rp->item[0];
gmp_printf(" %Zd\n", cp->data);
}
*/
        cell_ptr cpi, cpj;
        relmcount--;
        i=0; j=1;
        while (i<relmcount && j<rp->count - 1) {
          cpi = relm[i];
          cpj = rp->item[j];
          if (cpi->ind == cpj->ind) {
            mpz_mul(z1, z0, cpj->data);
            mpz_mod(z1, z1, order);
            int res = mpz_cmp(z1, cpi->data);
            if (!res) {
              memmove(&relm[i], &relm[i + 1], (relmcount - i - 1) * sizeof(cell_ptr));
              relm[relmcount - 1] = cpi;
              relmcount--;
              j++;
            } else if (res > 0) {
              mpz_sub(z1, order, z1);
              mpz_add(cpi->data, cpi->data, z1);
              i++;
              j++;
            } else {
              mpz_sub(cpi->data, cpi->data, z1);
              i++;
              j++;
            }
          } else if (cpi->ind > cpj->ind) {
            cpi = relm[relmcount];
            memmove(&relm[i + 1], &relm[i], (relmcount - i) * sizeof(cell_ptr));
            relm[i] = cpi;
            relmcount++;

            cpi->ind = cpj->ind;
            mpz_mul(z1, z0, cpj->data);
            mpz_mod(z1, z1, order);
            mpz_sub(cpi->data, order, z1);
            //cpi->data = order - ((u0 * cpj->data) % order);
            i++;
            j++;
          } else {
            i++;
          }
        }

        if (i == relmcount) {
          while (j < rp->count - 1) {
            cpi = relm[relmcount];
            cpj = rp->item[j];
            cpi->ind = cpj->ind;
            mpz_mul(z1, z0, cpj->data);
            mpz_mod(z1, z1, order);
            mpz_sub(cpi->data, order, z1);
            //cpi->data = order - ((u0 * cpj->data) % order);
            relmcount++;
            j++;
          }
        }

        cpj = rp->item[0];
        mpz_mul(z1, z0, cpj->data);
        mpz_mod(z1, z1, order);
        //u1 = (u0 * cpj->data) % order;
        if (mpz_cmp(km, z1) >= 0) {
          mpz_sub(km, km, z1);
        } else {
          mpz_sub(z1, order, z1);
          mpz_add(km, km, z1);
        }

/*
printf("after:");
for (i=0; i<relmcount; i++) {
  gmp_printf(" %u:%Zd", relm[i]->ind, relm[i]->data);
}
gmp_printf(" %Zd\n", km);
*/
      }
    }

  } while (count < r * fac->count);

  for (faci=0; faci<fac->count; faci++) {
    darray_t *row = matrix[faci];
    mpz_ptr order = fac->item[faci];
    for (i=1; i<r; i++) {
      darray_ptr rp = row[i];
      cell_ptr c0 = rp->item[0];
      for (j=1; j<rp->count-1; j++) {
        cell_ptr cp = rp->item[j];
        darray_ptr r2 = row[cp->ind];
        cell_ptr c2 = r2->item[0];
        mpz_mul(z0, cp->data, c2->data);
        mpz_sub(c0->data, c0->data, z0);
        mpz_mod(c0->data, c0->data, order);
      }
    }
  }

  mpz_ptr *tmp = pbc_malloc(sizeof(mpz_ptr) * fac->count);
  for (i=0; i<fac->count; i++) {
    tmp[i] = pbc_malloc(sizeof(mpz_t));
    mpz_init(tmp[i]);
    mpz_pow_ui(fac->item[i], fac->item[i], (unsigned int) mul->item[i]);
  }

  for (i=0; i<r; i++) {
    for (faci=0; faci<fac->count; faci++) {
      darray_t *row = matrix[faci];
      cell_ptr cp = row[i]->item[0];
      mpz_set(tmp[faci], cp->data);
    }
    CRT(ind[i], tmp, (mpz_ptr *) fac->item, fac->count);
  }

  for (i=0; i<fac->count; i++) {
    mpz_clear(tmp[i]);
  }
  pbc_free(tmp);

  for (faci=0; i<fac->count; faci++) {
    //similarly ''row'' refers to a list of rows
    darray_t *row = matrix[faci];
    for (j=0; j<r; j++) {
      darray_forall(row[j], delcell);
      darray_clear(row[j]);
    }
    pbc_free(matrix[faci]);
  }

  for (i=0; i<r; i++) {
    delcell(rel[i]);
    delcell(relm[i]);
  }

  pbc_free(prime);
  pbc_free(rel);
  pbc_free(relm);
  mpz_clear(k);
  mpz_clear(km);
  mpz_clear(z);
  mpz_clear(z0);
  mpz_clear(z1);
}

// Brute-force: does not use the fact that matrices are sparse.
void slow_index_calculus_step1(mpz_t *ind, int r, mpz_t g, mpz_t q,
    darray_ptr fac, darray_ptr mul) {
  int count = 0;
  int i, j;
  mpz_t z, z0, z1;
  //mpz_t rel[r + 1];
  //mpz_t relm[r + 1];
  mpz_t *rel = pbc_malloc(sizeof(mpz_t) * (r + 1));
  mpz_t *relm = pbc_malloc(sizeof(mpz_t) * (r + 1));
  unsigned int *prime = pbc_malloc(sizeof(unsigned int) * r);
  darray_t matrix;
  int faci;
  mpz_t k;
  int minfound[fac->count];

  for (i=0; i<r+1; i++) mpz_init(rel[i]);
  for (i=0; i<r+1; i++) mpz_init(relm[i]);

  mpz_init(k);
  mpz_init(z);
  mpz_init(z1);
  mpz_init(z0);

  darray_init(matrix);

  for (i=0; i<fac->count; i++) {
    darray_append(matrix, pbc_malloc(r * sizeof(mpz_t *)));
    minfound[i] = 0;
  }

  for (j=0; j<fac->count; j++) {
    mpz_t **row = matrix->item[j];
    for (i=0; i<r; i++) row[i] = NULL;
  }

  printf("building prime table...\n");
  prime[0] = 2;
  mpz_set_ui(z, 2);
  for (i=1; i<r; i++) {
    mpz_nextprime(z, z);
    prime[i] = mpz_get_ui(z);
  }
  printf("searching for r-smooth numbers\n");

  mpz_set_ui(z, 1);
  mpz_init(k);
  int try = 0;
  do {
    mpz_mul(z, z, g);
    mpz_mod(z, z, q);

    mpz_add_ui(k, k, 1);
    /*
    pbc_mpz_random(k, q);
    mpz_powm(z, g, k, q);
    */

    try++;

    mpz_set(z1, z);
    for (i=0; i<r; i++) {
      mpz_set_ui(rel[i], 0);
      while (mpz_divisible_ui_p(z1, prime[i])) {
        mpz_add_ui(rel[i], rel[i], 1);
        mpz_divexact_ui(z1, z1, prime[i]);
      }
    }
    if (mpz_cmp_ui(z1, 1)) {
      continue;
    }
    mpz_set(rel[r], k);

/*
    printf("found r-smooth number after %d tries\n", try);
    gmp_printf("g^%Zd = %Zd:", rel[r], z);
    for (i=0; i<r; i++) {
      if (mpz_sgn(rel[i])) {
        gmp_printf(" %u:%Zd", i, rel[i]);
      }
    }
    printf("\n");
*/

    try = 0;

    for (faci=0; faci<fac->count; faci++) {
      mpz_t **row = matrix->item[faci];
      mpz_ptr order = fac->item[faci];
      //gmp_printf("mod %Zd\n", order);
      for (i=0; i<r+1; i++) {
        mpz_mod(relm[i], rel[i], order);
      }

      for (;;) {
        /*
        for (i=0; i<r && !mpz_sgn(relm[i]); i++);
        if (i == r) {
          //printf("redundant relation\n");
          break;
        }
        */
        for (i=r-1; i>=0 && !mpz_sgn(relm[i]); i--);
        if (i < 0) {
          //printf("redundant relation\n");
          break;
        }
        if (i < minfound[faci]) {
          break;
        }
        mpz_set(z0, relm[i]);
        if (!row[i]) {
          row[i] = pbc_malloc(sizeof(mpz_t) * (r + 1));
          mpz_invert(z1, z0, order);
          for (j=0; j<r+1; j++) {
            mpz_init(row[i][j]);
            mpz_mul(row[i][j], z1, relm[j]);
            mpz_mod(row[i][j], row[i][j], order);
          }
          count++;
printf("%d / %d\n", count, r * fac->count);
/*
for (j=0; j<r; j++) {
  if (mpz_sgn(row[i][j])) {
    gmp_printf(" %d:%Zd", j, row[i][j]);
  }
}
gmp_printf(" %Zd\n", row[i][j]);
*/

          if (i == minfound[faci]) {
            do {
              if (!minfound[faci]) {
              printf("found log p_%d\n", minfound[faci]);
              gmp_printf("km = %Zd mod %Zd\n", row[i][r], order);
              }
              minfound[faci]++;
              if (minfound[faci] >= r) break;
            } while (row[minfound[faci]]);
          }
          break;
        }

        /*
        printf("before:");
        for (j=0; j<r; j++) {
          if (mpz_sgn(relm[j])) {
            gmp_printf(" %d:%Zd", j, relm[j]);
          }
        }
        gmp_printf(" %Zd\n", relm[j]);

        printf("sub %d:", i);
        for (j=0; j<r; j++) {
          if (mpz_sgn(row[i][j])) {
            gmp_printf(" %d:%Zd", j, row[i][j]);
          }
        }
        gmp_printf(" %Zd\n", row[i][j]);
        */

        for (j=0; j<r+1; j++) {
          mpz_mul(z1, z0, row[i][j]);
          mpz_sub(relm[j], relm[j], z1);
          mpz_mod(relm[j], relm[j], order);
        }

        /*
        printf("after:");
        for (j=0; j<r; j++) {
          if (mpz_sgn(relm[j])) {
            gmp_printf(" %d:%Zd", j, relm[j]);
          }
        }
        gmp_printf(" %Zd\n", relm[j]);
        */
      }
    }

  } while (count < r * fac->count);

  for (faci=0; faci<fac->count; faci++) {
    mpz_t **row = matrix->item[faci];
    mpz_ptr order = fac->item[faci];
    /*
    gmp_printf("mod %Zd:\n", order);
    for (i=0; i<r; i++) {
      for (j=0; j<r+1; j++) {
        gmp_printf(" %Zd", row[i][j]);
      }
      printf("\n");
    }
    printf("\n");
    */

    for (i=1; i<r; i++) {
      for (j=0; j<i; j++) {
        if (mpz_sgn(row[i][j])) {
          mpz_mul(z0, row[i][j], row[j][r]);
          mpz_sub(row[i][r], row[i][r], z0);
          mpz_mod(row[i][r], row[i][r], order);
        }
      }
    }
    /*
    for (i=r-2; i>=0; i--) {
      for (j=i+1; j<r; j++) {
        if (mpz_sgn(row[i][j])) {
          mpz_mul(z0, row[i][j], row[j][r]);
          mpz_sub(row[i][r], row[i][r], z0);
          mpz_mod(row[i][r], row[i][r], order);
        }
      }
    }
    */

    /*
    for (i=0; i<r; i++) {
      mpz_set(rel[i], row[i][r]);
      gmp_printf(" %Zd", row[i][r]);
      printf("\n");
    }
    */
  }

  mpz_ptr *tmp = pbc_malloc(sizeof(mpz_ptr) * fac->count);
  for (i=0; i<fac->count; i++) {
    tmp[i] = pbc_malloc(sizeof(mpz_t));
    mpz_init(tmp[i]);
    mpz_pow_ui(fac->item[i], fac->item[i], (unsigned int) mul->item[i]);
  }

  for (i=0; i<r; i++) {
    for (faci=0; faci<fac->count; faci++) {
      mpz_t **row = matrix->item[faci];
      mpz_set(tmp[faci], row[i][r]);
    }
    CRT(ind[i], tmp, (mpz_ptr *) fac->item, fac->count);
  }

  for (i=0; i<fac->count; i++) {
    mpz_clear(tmp[i]);
  }
  pbc_free(tmp);

  for (faci=0; faci<matrix->count; faci++) {
    mpz_t **row = matrix->item[faci];
    for (j=0; j<r; j++) {
      for (i=0; i<r+1; i++) {
        mpz_clear(row[j][i]);
      }
      pbc_free(row[j]);
    }
    pbc_free(row);
  }
  darray_clear(matrix);
  for (i=0; i<r+1; i++) mpz_clear(rel[i]);
  for (i=0; i<r+1; i++) mpz_clear(relm[i]);
  pbc_free(prime);
  pbc_free(rel);
  pbc_free(relm);
  mpz_clear(k);
  mpz_clear(z);
  mpz_clear(z0);
  mpz_clear(z1);

  printf("step 1 completed\n");
  for (i=0; i<r; i++) element_printf(" %Zd", ind[i]);
  printf("\n");
}

static void index_calculus_step2(mpz_t x, mpz_t *ind, int r,
    mpz_t g, mpz_t h, mpz_t q) {
  mpz_t prime;
  mpz_t s;
  mpz_t z, z1;
  mpz_t rel[r];
  int i;

  mpz_init(z);
  mpz_init(z1);
  mpz_init(s);
  mpz_init(prime);
  for (i=0; i<r; i++) mpz_init(rel[i]);

  mpz_set(z, h);

  for (;;) {
    mpz_mul(z, z, g);
    mpz_mod(z, z, q);
    mpz_add_ui(s, s, 1);

    mpz_set(z1, z);
    mpz_set_ui(prime, 1);
    for (i=0; i<r; i++) {
      mpz_set_ui(rel[i], 0);
      mpz_nextprime(prime, prime);
      while (mpz_divisible_p(z1, prime)) {
        mpz_add_ui(rel[i], rel[i], 1);
        mpz_divexact(z1, z1, prime);
      }
    }
    if (mpz_cmp_ui(z1, 1)) continue;
    element_printf("found r-smooth number on try #%Zd\n", s);
    mpz_set_ui(x, 0);
    for (i=0; i<r; i++) {
      mpz_mul(z, rel[i], ind[i]);
      mpz_add(x, x, z);
    }
    mpz_sub(x, x, s);
    mpz_sub_ui(z, q, 1);
    mpz_mod(x, x, z);
    break;
  }
}

static void mpzclear(void *p) {
  mpz_clear(p);
  pbc_free(p);
}

struct addfm_scope_var {
  darray_ptr fac, mul;
};

static int addfm(mpz_t f, unsigned int m, struct addfm_scope_var *v) {
  darray_append(v->fac, f);
  darray_append(v->mul, int_to_voidp(m));
  return 0;
}

void pbc_mpz_index_calculus(mpz_t x, mpz_t g, mpz_t h, mpz_t q) {
  int i, r;
  mpz_t q1, z0;

  mpz_init(q1);
  mpz_init(z0);

  mpz_sub_ui(q1, q, 1);
  mpz_setbit(z0, 6);

  darray_t fac, mul;
  darray_init(fac);
  darray_init(mul);
  struct addfm_scope_var v = {.fac = fac, .mul = mul};
  pbc_trial_divide((int(*)(mpz_t,unsigned,void*))addfm, &v, q1, z0);

  for (i=0; i<mul->count; i++) {
    unsigned int m = (unsigned int) mul->item[i];
    if (m != 1) {
      //TODO
      printf("p-adics not implemented yet\n");
      return;
    }
  }

  {
    double dq = mpz_get_d(q);
    //r = exp(sqrt(log(dq)*log(log(dq))));
    //printf("r = %d\n", r);
    r = exp(1.2 * sqrt(log(dq)));
    printf("r = %d\n", r);
  }
  mpz_t *ind = pbc_malloc(sizeof(mpz_t) * r);
  for (i=0; i<r; i++) mpz_init(ind[i]);

  if (is_gen(g, q, fac, mul)) {

    index_calculus_step1(ind, r, g, q, fac, mul);

    index_calculus_step2(x, ind, r, g, h, q);
  } else {
    mpz_t y, z;
    mpz_t d;

    mpz_init(d);
    mpz_init(y);
    mpz_init(z);
    do {
      pbc_mpz_random(z, q);
    } while (!is_gen(z, q, fac, mul));

    gmp_printf("new gen: %Zd\n", z);

    index_calculus_step1(ind, r, z, q, fac, mul);
    //slow_index_calculus_step1(ind, r, z, q, fac, mul);

    index_calculus_step2(x, ind, r, z, g, q);
    index_calculus_step2(y, ind, r, z, h, q);
    //want y / x mod q-1
    mpz_gcd(d, x, q1);
    mpz_divexact(q1, q1, d);
    mpz_divexact(x, x, d);
    //if y not divisible by d there is no solution
    mpz_divexact(y, y, d);
    mpz_invert(x, x, q1);
    mpz_mul(x, y, x);
    mpz_mod(x, x, q1);

    do {
      mpz_powm(z0, g, x, q);
      if (!mpz_cmp(z0, h)) {
        break;
      }
      mpz_add(x, x, q1);
      mpz_sub_ui(d, d, 1);
    } while (mpz_sgn(d));

    mpz_clear(d);
    mpz_clear(y);
    mpz_clear(z);
  }

  for (i=0; i<r; i++) mpz_clear(ind[i]);
  pbc_free(ind);

  darray_forall(fac, mpzclear);
  darray_clear(mul);
  darray_clear(fac);
  mpz_clear(q1);
  mpz_clear(z0);
}
