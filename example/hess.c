/*
  Hess ID-based signature.
  Based on papers "F. Hess. Efficient Identity Based Signature Schemes Based on Pairings. SAC 2002, LNCS 2595, Springer-Verlag, 2000"
  Contributed by Dmitry Kosolapov.
*/

#include <pbc.h>
#include <pbc_test.h>

int main(int argc, char **argv) {
  pairing_t pairing;
  double time1, time2;
  pbc_demo_pairing_init(pairing, argc, argv);

  element_t Qid, P, P1, Ppub, s, k, Did, r, v, u, t1, t3, t4, t5, t6, t7, t8;
  mpz_t t2;

  if (!pairing_is_symmetric(pairing)) pbc_die("pairing must be symmetric");

  mpz_init(t2);
  element_init_G1(P, pairing);
  element_init_G1(P1, pairing);
  element_init_G1(Qid, pairing);
  element_init_G1(Did, pairing);
  element_init_G1(Ppub, pairing);
  element_init_G1(t4, pairing);
  element_init_G1(t5, pairing);
  element_init_G1(u, pairing);

  element_init_Zr(s, pairing);
  element_init_Zr(k, pairing);
  element_init_Zr(v, pairing);
  element_init_Zr(t3, pairing);
  element_init_Zr(t8, pairing);

  element_init_GT(r, pairing);
  element_init_GT(t1, pairing);
  element_init_GT(t6, pairing);
  element_init_GT(t7, pairing);

  time1 = pbc_get_time();
  printf("Hess ID-based signature protocol\n");
  printf("KEYGEN\n");
  element_random(P);
  element_random(s);
  element_random(Qid);
  element_mul_zn(Ppub, P, s);
  element_mul_zn(Did, Qid, s);
  element_printf("Qid = %B\n", Qid);
  element_printf("P = %B\n", P);
  element_printf("Ppub = %B\n", Ppub);

  printf("SIGN\n");
  element_random(P1);
  element_random(k);
  element_pairing(t1, P1, P);
  element_pow_zn(r, t1, k);
  element_to_mpz(t2, r);

  //h3=h(m)*mpz(r);
  element_from_hash(t3, "Message", 7);
  element_mul_mpz(v, t3, t2);
  element_mul_zn(t4, Did, v);
  element_mul_zn(t5, P1, k);
  element_add(u, t4, t5);
  printf("Signature of message \"Message\" is:\n");
  element_printf("u = %B\n", u);
  element_printf("v = %B\n", v);

  printf("VERIFY\n");
  element_pairing(t6, u, P);
  element_neg(Ppub, Ppub);
  element_pairing(t7, Qid, Ppub);
  element_pow_zn(t7, t7, v);
  element_mul(r, t6, t7);
  element_to_mpz(t2, r);
  element_from_hash(t3, "Message", 7);
  element_mul_mpz(t8, t3, t2);
  element_printf("h3(m,r) = %B\n", t8);
  element_printf("v = %B\n", v);
  if (!element_cmp(t8, v)) {
    printf("Signature is valid!\n");
  } else {
    printf("Signature is invalid!\n");
  }
  time2 = pbc_get_time();
  printf("All time = %fs\n", time2 - time1);

  element_clear(P);
  element_clear(P1);
  element_clear(Qid);
  element_clear(Did);
  element_clear(Ppub);
  element_clear(t4);
  element_clear(t5);
  element_clear(u);
  element_clear(s);
  element_clear(k);
  element_clear(v);
  element_clear(t3);
  element_clear(t8);
  element_clear(r);
  element_clear(t1);
  element_clear(t6);
  element_clear(t7);
  pairing_clear(pairing);

  return 0;
}
