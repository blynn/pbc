// Mutliexponentiation benchmark and test.

#include <string.h>
#include "pbc.h"
#include "pbc_test.h"

int main(int argc, char **argv) {
  pairing_t pairing;
  element_t g1, u1, up1, g2, u2, up2, r;
  mpz_t r_mpz;
  element_pp_t g1_pp, g2_pp;
  double t0, t1;
  int i, n;

  printf("reading pairing from stdin...\n");
  pbc_demo_pairing_init(pairing, argc, argv);

  element_init(r, pairing->Zr);
  element_init(g1, pairing->G1);
  element_init(u1, pairing->G1);
  element_init(up1, pairing->G1);
  element_init(g2, pairing->G2);
  element_init(u2, pairing->G2);
  element_init(up2, pairing->G2);

  element_random(r);
  element_random(g1);
  element_random(g2);

  mpz_init(r_mpz);
  element_to_mpz(r_mpz, r);

  element_pp_init(g1_pp, g1);
  element_pp_init(g2_pp, g2);

  n = 100;
  t0 = pbc_get_time();
  for (i=0; i<n; i++) {
    element_pow_mpz(u1, g1, r_mpz);
  }
  t1 = pbc_get_time();
  printf("G1 exp:\t\t%fs\n", t1 - t0);

  n = 100;
  t0 = pbc_get_time();
  for (i=0; i<n; i++) {
    element_pow_mpz(u2, g2, r_mpz);
  }
  t1 = pbc_get_time();
  printf("G2 exp:\t\t%fs\n", t1 - t0);

  n = 100;
  t0 = pbc_get_time();
  for (i=0; i<n; i++) {
    element_pp_pow(up1, r_mpz, g1_pp);
  }
  t1 = pbc_get_time();
  printf("G1 pp exp:\t%fs\n", t1 - t0);

  n = 100;
  t0 = pbc_get_time();
  for (i=0; i<n; i++) {
    element_pp_pow(up2, r_mpz, g2_pp);
  }
  t1 = pbc_get_time();
  printf("G2 pp exp:\t%fs\n", t1 - t0);

  if (element_cmp(u1, up1)) {
    printf("Oops 1!\n");
  }
  if (element_cmp(u2, up2)) {
    printf("Oops 2!\n");
  }

  mpz_clear(r_mpz);
  element_clear(g1);
  element_clear(u1);
  element_clear(up1);
  element_clear(g2);
  element_clear(u2);
  element_clear(up2);
  element_clear(r);
  element_pp_clear(g1_pp);
  element_pp_clear(g2_pp);
  pairing_clear(pairing);

  return 0;
}
