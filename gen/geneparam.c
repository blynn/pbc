// Generate type E pairings.
// Usage:
//   geneparam [RBITS [QBITS]]
//
// RBITS
//   The number of bits in r, the order of the subgroup G1. Default is 160.
// QBITS
//   The number of bits in q, the order of the full group. Default is 1024.

#include "pbc.h"

int main(int argc, char **argv) {
  int rbits = argc > 1 ? atoi(argv[1]) : 160;
  int qbits = argc > 2 ? atoi(argv[2]) : 1024;

  pbc_param_t par;
  pbc_param_init_e_gen(par, rbits, qbits);
  pbc_param_out_str(stdout, par);
  pbc_param_clear(par);
  return 0;
}
