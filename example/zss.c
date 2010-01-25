/*
  ZSS Short Signature Scheme from Bilinear Pairing.
  Based on papers "F. Zhang, R. Safavi-Naini and W. Susilo. An Efficient Signature Scheme from Bilinear Pairings and it's Applications. PKC 2004".
  Contributed by Dmitry Kosolapov.
*/

#include <pbc.h>
#include <pbc_test.h>

int main(int argc, char **argv) {
  pairing_t pairing;
  pbc_demo_pairing_init(pairing, argc, argv);
  if (!pairing_is_symmetric(pairing)) pbc_die("pairing must be symmetric");
  double time1, time2;
  element_t P, Ppub, x, S, H, t1, t2, t3, t4;
  element_init_Zr(x, pairing);
  element_init_Zr(H, pairing);
  element_init_Zr(t1, pairing);

  element_init_G1(S, pairing);
  element_init_G1(P, pairing);
  element_init_G1(Ppub, pairing);
  element_init_G1(t2, pairing);

  element_init_GT(t3, pairing);
  element_init_GT(t4, pairing);

  printf("ZSS short signature schema\n");
  printf("KEYGEN\n");
  time1 = pbc_get_time();
  element_random(x);
  element_random(P);
  element_mul_zn(Ppub, P, x);
  element_printf("P = %B\n", P);
  element_printf("x = %B\n", x);
  element_printf("Ppub = %B\n", Ppub);

  printf("SIGN\n");
  element_from_hash(H, "Message", 7);
  element_add(t1, H, x);
  element_invert(t1, t1);
  element_mul_zn(S, P, t1);
  printf("Signature of message \"Message\" is:\n");
  element_printf("S = %B\n", S);

  printf("VERIFY\n");
  element_from_hash(H, "Message", 7);
  element_mul_zn(t2, P, H);
  element_add(t2, t2, Ppub);
  element_pairing(t3, t2, S);
  element_pairing(t4, P, P);
  element_printf("e(H(m)P + Ppub, S) = %B\n", t3);
  element_printf("e(P, P) = %B\n", t4);
  if (!element_cmp(t3, t4)) printf("Signature is valid\n");
  else printf("Signature is invalid\n");
  time2 = pbc_get_time();
  printf("All time = %fs\n", time2 - time1);
  element_clear(P);
  element_clear(Ppub);
  element_clear(x);
  element_clear(S);
  element_clear(H);
  element_clear(t1);
  element_clear(t2);
  element_clear(t3);
  element_clear(t4);
  pairing_clear(pairing);

  return 0;
}
