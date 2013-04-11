// Generate MNT curve(s) for a given D.

#include <stdint.h> // for intptr_t
#include "pbc.h"
#include "pbc_utils.h"  // for UNUSED_VAR

int generate(pbc_cm_t cm, void *data) {
  UNUSED_VAR(data);
  pbc_param_t param;
  pbc_info("gendparam: computing Hilbert polynomial and finding roots...");
  pbc_param_init_d_gen(param, cm);
  pbc_info("gendparam: bits in q = %zu\n", mpz_sizeinbase(cm->q, 2));
  pbc_param_out_str(stdout, param);
  pbc_param_clear(param);
  return 1;
}

int main(int argc, char **argv) {
  int D = 9563;

  if (argc > 1) {
    int m;
    D = atoi(argv[1]);
    m = D % 4;
    if (D <= 0 || m == 1 || m == 2) {
      pbc_die("D must be 0 or 3 mod 4 and positive");
    }
  }
  pbc_info("Using D = %d\n", D);

  if (!pbc_cm_search_d(generate, NULL, D, 500)) {
    pbc_die("no suitable curves for this D");
  }
  return 0;
}
