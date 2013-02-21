// List discriminant and bits in r and q for type G pairings that may be
// suitable for cryptography.

#include "pbc.h"

int show(pbc_cm_t cm, void *data) {
  unsigned int D = * (unsigned *) data;
  int qbits, rbits;
  qbits = mpz_sizeinbase(cm->q, 2);
  rbits = mpz_sizeinbase(cm->r, 2);
  printf("%d, %d, %d\n", D, qbits, rbits);
  fflush(stdout);
  return 0;
}

void try(int tryD) {
  pbc_cm_search_g(show, &tryD, tryD, 500);
}

int main(int argc, char **argv) {
  unsigned int D = 0;

  if (argc > 1) {
    D = atoi(argv[1]);
    if (D % 120) {
      pbc_die("D must be multiple of 120");
    }
  }

  printf("D < %u, bits in q, bits in r\n", 1000000000);
  while (D < 1000000000) {
    try(D + 43);
    try(D + 67);
    D+=120;
  }

  return 0;
}
