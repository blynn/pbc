// Prints Hilbert polynomials H_D(X)
//
// Usage: hilbertpoly [LOWER [UPPER]]
//
// LOWER:
//   Lower limit of D. Defaults to 3.
// UPPER:
//   Upper limit of D. Defaults to LOWER.
//
// e.g. $ hilbertpoly 3 1000000
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h> //for atoi, exit
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_field.h"
#include "pbc_darray.h"
#include "pbc_poly.h"
#include "pbc_hilbert.h"

int main(int argc, char **argv) {
  int n;
  int i;
  darray_t coefflist;
  void xpow(int degree) {
    if (degree == 1) {
      printf("X");
    } else if (degree) {
      printf("X^%d", degree);
    }
  }

  int D      = argc > 1 ? atoi(argv[1]) : 3;
  if (D <= 0) pbc_die("D must be positive.");

  int Dlimit = argc > 2 ? atoi(argv[2]) : D;

  for(; D <= Dlimit; D++) {
    int m = D % 4;
    if (m == 1 || m == 2) continue;
    printf("D = %d\n", D);

    darray_init(coefflist);
    pbc_hilbert(coefflist, D);

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
    pbc_hilbert_clear(coefflist);
    darray_clear(coefflist);
  }

  return 0;
}
