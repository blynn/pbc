#include <stdio.h>
#include <stdlib.h> //for atoi, exit
#include <gmp.h>
#include "field.h"
#include "darray.h"
#include "poly.h"
#include "hilbert.h"

int main(int argc, char **argv)
{
    int n;
    int i;
    int D = 3;
    darray_t coefflist;
    void xpow(int degree) {
	if (degree == 1) {
	    printf("X");
	} else if (degree) {
	    printf("X^%d", degree);
	}
    }

    if (argc > 1) {
	int m;
	D = atoi(argv[1]);
	m = D % 4;
	if (D <= 0 || m == 1 || m == 2) {
	    fprintf(stderr, "D must be 0 or 3 mod 4 and positive\n");
	    exit(1);
	}
    }
    printf("D = %d\n", D);

    darray_init(coefflist);

    hilbert_poly(coefflist, D);

    n = coefflist->count;
    printf(" ");
    xpow(n - 1);
    printf("\n");
    for (i=n-2; i>=0; i--) {
	if (mpz_sgn((mpz_ptr) coefflist->item[i]) >= 0) {
	    printf("+");
	}
	mpz_out_str(stdout, 0, coefflist->item[i]);
	xpow(i);
	printf("\n");
    }

    return 0;
}
