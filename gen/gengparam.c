// Generate Freeman curves with a given discriminant.
#include "pbc.h"

void generate(cm_info_t cm, void *data) {
  UNUSED_VAR(data);
  pbc_param_t param;

  pbc_info("gengparam: computing Hilbert polynomial and finding roots...");
  pbc_param_init_g_gen(param, cm);
  pbc_info("gengparam: bits in q = %zu", mpz_sizeinbase(cm->q, 2));
  pbc_param_out_str(stdout, param);
  pbc_param_clear(param);
}

int main(int argc, char **argv) {
  int D = 35707;

  if (argc > 1) {
    int m;
    D = atoi(argv[1]);
    m = D % 120;
    if (D <= 0 || (m != 43 && m != 67)) {
      pbc_die("D must be 43 or 67 mod 120 and positive");
    }
  }
  pbc_info("Using D = %d", D);

  if (!find_freeman_curve(generate, NULL, D, 500)) {
    pbc_die("No suitable curves for this D");
  }

  return 0;
}
