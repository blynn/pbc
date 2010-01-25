// Compares dedicated multipairing (product of pairings) routine with naive
// method.
#include <pbc.h>
#include "pbc_test.h"

int main(int argc, char **argv) {
  enum { K = 5 };
  pairing_t pairing;
  element_t x[K], y[K], r, r2, tmp;
  int i, n;
  double t0, t1, ttotal, ttotalm;

  pbc_demo_pairing_init(pairing, argc, argv);

  for(i = 0; i < K; i++) {
    element_init_G1(x[i], pairing);
    element_init_G2(y[i], pairing);
  }
  element_init_GT(r, pairing);
  element_init_GT(r2, pairing);
  element_init_GT(tmp, pairing);

  n = 10;
  ttotal = 0.0;
  ttotalm = 0.0;
  for (i=0; i<n; i++) {
    int j;
    for(j = 0; j < K; j++) {
      element_random(x[j]);
      element_random(y[j]);
    }

    t0 = pbc_get_time();
    element_prod_pairing(r, x, y, K);
    t1 = pbc_get_time();
    ttotalm += t1 - t0;

    t0 = pbc_get_time();
    element_pairing(r2, x[0], y[0]);
    for(j = 1; j < K; j++) {
      element_pairing(tmp, x[j], y[j]);
      element_mul(r2, r2, tmp);
    }
    t1 = pbc_get_time();
    ttotal += t1 - t0;

    element_printf("e(x,y) = %B\n", r);
    EXPECT(!element_cmp(r, r2));
  }
  printf("average pairing time = %f\n", ttotal / n);
  printf("average multi-pairing time = %f\n", ttotalm / n);

  for(i = 0; i < K; i++) {
    element_clear(x[i]);
    element_clear(y[i]);
  }
  element_clear(r);
  element_clear(r2);

  pairing_clear(pairing);
  return 0;
}
