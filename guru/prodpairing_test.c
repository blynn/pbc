// Check product of pairings works for F pairings when initialized via
// pairing_init_pbc_param().
//
// By Michael Adjedj, Ben Lynn.
#include "pbc.h"
#include "pbc_test.h"

int main(void) {
  pbc_param_t param;

  pbc_param_init_f_gen(param, 200);
  pairing_t pairing;
  pairing_init_pbc_param(pairing, param);

  element_t P[2], Q[2], res, tmp, tmp2;

  element_init_G1(P[0], pairing);  element_random(P[0]);
  element_init_G1(P[1], pairing);  element_random(P[1]);

  element_init_G2(Q[0], pairing);  element_random(Q[0]);
  element_init_G2(Q[1], pairing);  element_random(Q[1]);

  element_init_GT(res, pairing);
  element_init_GT(tmp, pairing);
  element_init_GT(tmp2, pairing);

  element_prod_pairing(res, P, Q, 2);
  element_pairing(tmp, P[0], Q[0]);
  element_pairing(tmp2, P[1], Q[1]);
  element_mul(tmp, tmp, tmp2);
  EXPECT(!element_cmp(res, tmp));

  element_clear(P[0]);
  element_clear(P[1]);
  element_clear(Q[0]);
  element_clear(Q[1]);
  element_clear(res);
  element_clear(tmp);
  element_clear(tmp2);

  pairing_clear(pairing);
  pbc_param_clear(param);
  return 0;
}
