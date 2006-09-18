//Find MNT curves with embedding degree 6

#include <stdlib.h>
#include "mnt.h"

void cm_info_init(cm_info_t cm)
{
    mpz_init(cm->q);
    mpz_init(cm->r);
    mpz_init(cm->h);
    mpz_init(cm->n);
}

void cm_info_clear(cm_info_t cm)
{
    mpz_clear(cm->q);
    mpz_clear(cm->r);
    mpz_clear(cm->h);
    mpz_clear(cm->n);
}

int mnt_step2(darray_ptr L, unsigned int D, mpz_t U)
{
    int d;
    mpz_t n, l, q;
    mpz_t p;
    mpz_t r, cofac;
    cm_info_ptr cm;

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
	return 1;
    }

    mpz_divexact_ui(l, l, 3);
    mpz_init(q);

    mpz_mul(q, l, l);
    mpz_add_ui(q, q, 1);
    if (!mpz_probab_prime_p(q, 10)) {
	mpz_clear(q);
	mpz_clear(l);
	return 1;
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
	    return 1;
	}
    }
}

    cm = malloc(sizeof(cm_info_t));
    cm_info_init(cm);
    cm->k = 6;
    cm->D = D;
    mpz_set(cm->q, q);
    mpz_set(cm->r, r);
    mpz_set(cm->h, cofac);
    mpz_set(cm->n, n);
    darray_append(L, cm);

    mpz_clear(cofac);
    mpz_clear(r);
    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(l);
    mpz_clear(n);
    return 0;
}

int find_mnt6_curve(darray_t L, unsigned int D, unsigned int bitlimit)
{
    //first solve DV^2 = 3l^2 - 2l + 3
    //use variable subst U = 3l - 1 to transform equation
    //U^2 - 3DV^2 = -8
    //so 3D cannot be a square
    //(the only squares that differ by 8 are 1 and 9,
    //which we get if U=V=1, D=3, but then l is not an integer)

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
    int found_count = 0;

    mpz_t D3;
    
    darray_t listp, listq;
    mpz_ptr zptr;

    mpz_init(a0); mpz_init(twice_a0); mpz_init(a1);
    mpz_init(P0); mpz_init(P1);
    mpz_init(Q0); mpz_init(Q1);
    mpz_init(p0); mpz_init(p1); mpz_init(pnext);
    mpz_init(q0); mpz_init(q1); mpz_init(qnext);
    mpz_init(D3);
    mpz_set_ui(D3, 3 * D);

    darray_init(listp);
    darray_init(listq);

    if (mpz_perfect_square_p(D3)) {
	return 0;
    }
    mpz_sqrt(a0, D3);
    mpz_set_ui(P0, 0);
    mpz_set_ui(Q0, 1);

    mpz_set(P1, a0);
    mpz_mul(Q1, a0, a0);
    mpz_sub(Q1, D3, Q1);
    mpz_add(a1, a0, P1);
    mpz_tdiv_q(a1, a1, Q1);

    mpz_add(twice_a0, a0, a0);

    mpz_set(p0, a0);
    mpz_set_ui(q0, 1);
    mpz_mul(p1, a0, a1);
    mpz_add_ui(p1, p1, 1);
    mpz_set(q1, a1);

    //mpz_out_str(stdout, 10, p0);
    //printf(" / ");
    //mpz_out_str(stdout, 10, q0);
    //printf("\n");

    d = -1;
    for(;;) {
	//printf("Q: ");
	if (d == -1) {
	    //printf("-");
	    if (!mpz_cmp_ui(Q1, 8)) {
		zptr = (mpz_ptr) malloc(sizeof(mpz_t));
		mpz_init(zptr);
		mpz_set(zptr, p0);
		darray_append(listp, zptr);
		zptr = (mpz_ptr) malloc(sizeof(mpz_t));
		mpz_init(zptr);
		mpz_set(zptr, q0);
		darray_append(listq, zptr);
	    } else if (!mpz_cmp_ui(Q1, 2)) {
		zptr = (mpz_ptr) malloc(sizeof(mpz_t));
		mpz_init(zptr);
		mpz_mul_ui(zptr, p0, 2);
		darray_append(listp, zptr);
		zptr = (mpz_ptr) malloc(sizeof(mpz_t));
		mpz_init(zptr);
		mpz_mul_ui(zptr, q0, 2);
		darray_append(listq, zptr);
	    }
	}
	//mpz_out_str(stdout, 10, Q1);
	//printf("\n");
	//mpz_out_str(stdout, 10, p1);
	//printf(" / ");
	//mpz_out_str(stdout, 10, q1);
	//printf("\n");
	if (!mpz_cmp(twice_a0, a1)) break;
	//compute more of the continued fraction expansion
	mpz_set(P0, P1);
	mpz_mul(P1, a1, Q1);
	mpz_sub(P1, P1, P0);
	mpz_set(Q0, Q1);
	mpz_mul(Q1, P1, P1);
	mpz_sub(Q1, D3, Q1);
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

    {
	int i, n;
	n = listp->count;
	if (n) for (;;) {
	    for (i=0; i<n; i++) {
		/*
		mpz_out_str(stdout, 0, listp->item[i]);
		printf(", ");
		mpz_out_str(stdout, 0, listq->item[i]);
		printf("\n");
		*/
		if (!mnt_step2(L, D, listp->item[i])) found_count++;
		//compute next solution as follows
		//if p, q is current solution
		//compute new solution p', q' via
		//(p + q sqrt{3D})(t + u sqrt{3D}) = p' + q' sqrt(3D)
		//where t, u is min. solution to Pell equation
		//can use a0, p1, q1 as temp variables now
		mpz_mul(p1, p0, listp->item[i]);
		mpz_mul(q1, q0, listq->item[i]);
		mpz_mul(q1, q1, D3);
		mpz_add(p1, p1, q1);
		if (2 * mpz_sizeinbase(p1, 2) > bitlimit + 10) goto toobig;
		mpz_mul(a0, p0, listq->item[i]);
		mpz_mul(q1, q0, listp->item[i]);
		mpz_add(a0, a0, q1);
		mpz_set(listp->item[i], p1);
		mpz_set(listq->item[i], a0);
	    }
	}
    }
toobig:

    mpz_clear(a0); mpz_clear(twice_a0); mpz_clear(a1);
    mpz_clear(P0); mpz_clear(P1);
    mpz_clear(Q0); mpz_clear(Q1);
    mpz_clear(p0); mpz_clear(p1); mpz_clear(pnext);
    mpz_clear(q0); mpz_clear(q1); mpz_clear(qnext);
    mpz_clear(D3);
    return found_count;
}
