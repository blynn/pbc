#include <pbc.h>
#include "pbc_test.h"

void time_pairing(pairing_ptr pairing) {
  int i, n;
  double t0, t1, ttotal, ttotalpp;
  pairing_pp_t pp;
  element_t x, y, r, r2;
  element_init_G1(x, pairing);
  element_init_G2(y, pairing);
  element_init_GT(r, pairing);
  element_init_GT(r2, pairing);

  n = 10;
  ttotal = 0.0;
  ttotalpp = 0.0;
  for (i=0; i<n; i++) {
    element_random(x);
    element_random(y);

    pairing_pp_init(pp, x, pairing);
    t0 = pbc_get_time();
    pairing_pp_apply(r, y, pp);
    t1 = pbc_get_time();
    ttotalpp += t1 - t0;
    pairing_pp_clear(pp);

    t0 = pbc_get_time();
    element_pairing(r2, x, y);
    t1 = pbc_get_time();
    ttotal += t1 - t0;

    //element_printf("x = %B\n", x);
    //element_printf("y = %B\n", y);
    //element_printf("e(x,y) = %B\n", r);
    if (element_cmp(r, r2)) {
      printf("BUG!\n");
      exit(1);
    }
  }
  printf("average pairing time = %f\n", ttotal / n);
  printf("average pairing time (preprocessed) = %f\n", ttotalpp / n);

  element_clear(x);
  element_clear(y);
  element_clear(r);
  element_clear(r2);
}

int main(int argc, char **argv) {
  pairing_t pairing;

  pbc_demo_pairing_init(pairing, argc, argv);

  printf("Miller's algorithm\n");
  time_pairing(pairing);

  pairing_option_set(pairing, "method", "shipsey-stange");
  printf("Shipsey-Stange algorithm\n");
  time_pairing(pairing);

  pairing_clear(pairing);

  return 0;
}
