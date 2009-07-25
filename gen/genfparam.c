// Generate type F pairings.
// Usage:
//   genaparam [BITS]
//
// BITS
//   The number of bits in r, the order of the subgroup G1. Default is 160.

#include "pbc.h"

int main(int argc, char **argv) {
  int bits = 160;
  if (argc > 1) {
    bits = atoi(argv[1]);
    if (bits < 1) {
      pbc_die("Usage: %s [BITS]", argv[0]);
    }
  }
  pbc_param_t fp;
  pbc_param_init_f_gen(fp, bits);
  pbc_param_out_str(stdout, fp);
  pbc_param_clear(fp);

  return 0;
}
