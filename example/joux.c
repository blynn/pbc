/*
  Joux one round protocol for tripartite Diffie-Hellman
  Based on papers "A. Joux. A One Round Protocol for Tripartie Diffie-Hellman. Proceedings of ANTS 4. LNCS 1838, pp. 385-394, 2000."
  Contributed by Dmitry Kosolapov.
*/

#include <pbc.h>
#include <pbc_test.h>

int main(int argc, char **argv) {
  pairing_t pairing;
  double time1, time2;
  element_t P, a, b, c, Ka, Kb, Kc, t1, t2, t3, t4, t5, t6;
  pbc_demo_pairing_init(pairing, argc, argv);
  if (!pairing_is_symmetric(pairing)) pbc_die("pairing must be symmetric");

  element_init_G1(P, pairing);
  element_init_G1(t1, pairing);
  element_init_G1(t2, pairing);
  element_init_G1(t3, pairing);

  element_init_Zr(a, pairing);
  element_init_Zr(b, pairing);
  element_init_Zr(c, pairing);

  element_init_GT(t4, pairing);
  element_init_GT(t5, pairing);
  element_init_GT(t6, pairing);
  element_init_GT(Ka, pairing);
  element_init_GT(Kb, pairing);
  element_init_GT(Kc, pairing);

  time1 = pbc_get_time();
  printf("Joux key agreement between A, B and C.\n");
  element_random(P);
  element_random(a);
  element_random(b);
  element_random(c);
  element_mul_zn(t1, P, a);
  printf("A sends B and C: aP\n");
  element_printf("aP = %B\n", t1);
  element_mul_zn(t2, P, b);
  printf("B sends A and C: bP\n");
  element_printf("bP = %B\n", t2);
  element_mul_zn(t3, P, c);
  printf("C sends A and B: cP\n");
  element_printf("cP = %B\n", t3);

  element_pairing(t4, t2, t3);
  element_pow_zn(Ka, t4, a);
  element_printf("Ka = %B\n", Ka);
  element_pairing(t5, t1, t3);
  element_pow_zn(Kb, t5, b);
  element_printf("Kb = %B\n", Kb);
  element_pairing(t6, t1, t2);
  element_pow_zn(Kc, t6, c);
  element_printf("Kc = %B\n", Kc);

  printf("Shared key K = Ka = Kb = Kc\n");
  time2 = pbc_get_time();
  printf("All time = %fs\n", time2 - time1);


  element_clear(P);
  element_clear(a);
  element_clear(b);
  element_clear(c);
  element_clear(Ka);
  element_clear(Kb);
  element_clear(Kc);
  element_clear(t1);
  element_clear(t2);
  element_clear(t3);
  element_clear(t4);
  element_clear(t5);
  element_clear(t6);
  pairing_clear(pairing);

  return 0;
}
