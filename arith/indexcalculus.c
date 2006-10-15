#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include "pbc.h"
#include "random.h"
#include "utils.h"

static int is_gen(mpz_t x, mpz_t q, darray_ptr fac, darray_ptr mul)
{
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

static void CRT(mpz_t x, mpz_ptr *v, mpz_ptr *m, int t)
//Garner's Algorithm
//see Algorithm 14.71, Handbook of Cryptography
{
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

static void index_calculus_step1(mpz_t *index, int r, mpz_t g, mpz_t q,
	darray_ptr fac, darray_ptr mul)
{
    int count = 0;
    int i, j;
    mpz_t prime;
    mpz_t z, z0, z1;
    //mpz_t *row[r];
    mpz_t rel[r + 1];
    mpz_t relm[r + 1];
    darray_t matrix;
    int rowi;
    mpz_t k;

    for (i=0; i<r+1; i++) mpz_init(rel[i]);
    for (i=0; i<r+1; i++) mpz_init(relm[i]);

    mpz_init(prime);
    mpz_init(k);
    mpz_init(z);
    mpz_init(z1);
    mpz_init(z0);

    darray_init(matrix);

    for (i=0; i<fac->count; i++) {
	darray_append(matrix, malloc(r * sizeof(mpz_t *)));
    }

    for (j=0; j<fac->count; j++) {
	mpz_t **row = matrix->item[j];
	for (i=0; i<r; i++) row[i] = NULL;
    }

    mpz_set_ui(z, 1);
    mpz_init(k);
    do {
	mpz_mul(z, z, g);
	mpz_mod(z, z, q);
	mpz_add_ui(k, k, 1);
	if (!mpz_cmp_ui(z, 1)) {
	    fprintf(stderr, "went through whole group!\n");
	    break;
	    exit(0);
	}
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
	mpz_set(rel[r], k);
	/*
	printf("found r-smooth number\n");

	gmp_printf("z = %Zd: ", z);
	for (i=0; i<r+1; i++) {
	    gmp_printf(" %Zd", rel[i]);
	}
	printf("\n");
	*/

	for (rowi=0; rowi<fac->count; rowi++) {
	    mpz_t **row = matrix->item[rowi];
	    mpz_ptr order = fac->item[rowi];
	    //gmp_printf("mod %Zd\n", order);
	    for (i=0; i<r+1; i++) {
		mpz_mod(relm[i], rel[i], order);
	    }

	    for (;;) {
		for (i=0; i<r && !mpz_sgn(relm[i]); i++);
		if (i == r) {
		    //printf("redundant relation\n");
		    break;
		}
		mpz_set(z0, relm[i]);
		if (!row[i]) {
		    row[i] = malloc(sizeof(mpz_t) * (r + 1));
		    mpz_invert(z1, z0, order);
		    for (j=0; j<r+1; j++) {
			mpz_init(row[i][j]);
			mpz_mul(row[i][j], z1, relm[j]);
			mpz_mod(row[i][j], row[i][j], order);
		    }
		    count++;
printf("%d / %d\n", count, r * fac->count);
		    break;
		}
		/*
		printf("before:");
		for (j=0; j<r+1; j++) {
		    gmp_printf(" %Zd", relm[j]);
		}
		printf("\n");
		*/

		for (j=0; j<r+1; j++) {
		    mpz_mul(z1, z0, row[i][j]);
		    mpz_sub(relm[j], relm[j], z1);
		    mpz_mod(relm[j], relm[j], order);
		}

		/*
		printf("after:");
		for (j=0; j<r+1; j++) {
		    gmp_printf(" %Zd", relm[j]);
		}
		printf("\n");
		*/
	    }
	}

    } while (count < r * fac->count);

    for (rowi=0; rowi<fac->count; rowi++) {
	mpz_t **row = matrix->item[rowi];
	mpz_ptr order = fac->item[rowi];
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

	for (i=r-2; i>=0; i--) {
	    for (j=i+1; j<r; j++) {
		if (mpz_sgn(row[i][j])) {
		    mpz_mul(z0, row[i][j], row[j][r]);
		    mpz_sub(row[i][r], row[i][r], z0);
		    mpz_mod(row[i][r], row[i][r], order);
		}
	    }
	}

	/*
	for (i=0; i<r; i++) {
	    mpz_set(rel[i], row[i][r]);
	    gmp_printf(" %Zd", row[i][r]);
	    printf("\n");
	}
	*/
    }

    mpz_ptr *tmp = malloc(sizeof(mpz_ptr) * fac->count);
    for (i=0; i<fac->count; i++) {
	tmp[i] = malloc(sizeof(mpz_t));
	mpz_init(tmp[i]);
	mpz_pow_ui(fac->item[i], fac->item[i], (unsigned int) mul->item[i]);
    }

    for (i=0; i<r; i++) {
	for (rowi=0; rowi<fac->count; rowi++) {
	    mpz_t **row = matrix->item[rowi];
	    mpz_set(tmp[rowi], row[i][r]);
	}
	CRT(index[i], tmp, (mpz_ptr *) fac->item, fac->count);
    }

    for (i=0; i<fac->count; i++) {
	mpz_clear(tmp[i]);
    }
    free(tmp);

    for (rowi=0; rowi<matrix->count; rowi++) {
	mpz_t **row = matrix->item[rowi];
	for (j=0; j<r; j++) {
	    for (i=0; i<r+1; i++) {
		mpz_clear(row[j][i]);
	    }
	    free(row[j]);
	}
	free(row);
    }
    darray_clear(matrix);
    for (i=0; i<r+1; i++) mpz_clear(rel[i]);
    for (i=0; i<r+1; i++) mpz_clear(relm[i]);
    mpz_clear(k);
    mpz_clear(z);
    mpz_clear(z0);
    mpz_clear(z1);
    mpz_clear(prime);
}

static void index_calculus_step2(mpz_t x, mpz_t *index, int r,
	mpz_t g, mpz_t h, mpz_t q)
{
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
	    mpz_mul(z, rel[i], index[i]);
	    mpz_add(x, x, z);
	}
	mpz_sub(x, x, s);
	mpz_sub_ui(z, q, 1);
	mpz_mod(x, x, z);
	break;
    }
}

void index_calculus_dlog(mpz_t x, mpz_t g, mpz_t h, mpz_t q)
{
    darray_t fac, mul;
    int i, r;
    mpz_t q1, z0;

    void mpzclear(void *p)
    {
	mpz_clear(p);
	free(p);
    }

    darray_init(fac);
    darray_init(mul);

    mpz_init(q1);
    mpz_init(z0);

    mpz_sub_ui(q1, q, 1);
    mpz_setbit(z0, 6);
    trial_divide(fac, mul, q1, z0);

    for (i=0; i<mul->count; i++) {
	unsigned int m = (unsigned int) mul->item[i];
	if (m != 1) {
	    //TODO
	    printf("p-adics not implemented yet\n");
	    return;
	}
    }

    //r = 1000;
    {
	double dq = mpz_get_d(q);
	//r = exp(2 * sqrt(log(dq)*log(log(dq))));
	r = exp(sqrt(log(dq)));
	printf("r = %d\n", r);
    }
    mpz_t *index = malloc(sizeof(mpz_t) * r);
    for (i=0; i<r; i++) mpz_init(index[i]);

    if (is_gen(g, q, fac, mul)) {

	index_calculus_step1(index, r, g, q, fac, mul);

	printf("step 1 completed\n");
	//for (i=0; i<r; i++) element_printf(" %Zd", index[i]);
	//printf("\n");

	index_calculus_step2(x, index, r, g, h, q);
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

	index_calculus_step1(index, r, z, q, fac, mul);
	index_calculus_step2(x, index, r, z, g, q);
	index_calculus_step2(y, index, r, z, h, q);
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

    for (i=0; i<r; i++) mpz_clear(index[i]);
    free(index);

    darray_forall(fac, mpzclear);
    darray_clear(mul);
    darray_clear(fac);
    mpz_clear(q1);
    mpz_clear(z0);
}
