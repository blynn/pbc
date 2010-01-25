/*
  Paterson ID-based signature.
  Based on papers "K. G. Paterson. ID-Based Signatures from Pairings on Elliptic Curvers. Electron. Lett., Vol. 38". Available at http://eprint.iacr.org/2002/004."
  Contributed by Dmitry Kosolapov.
*/

#include <pbc.h>
#include <pbc_test.h>

int main(int argc, char **argv) {
  pairing_t pairing;
  double time1, time2;
  element_t Ppub, s, P, R, k, S, Did, Qid, t1, t2, t4, t5, t6, t7, t8,
    t9, t10, t11;
  mpz_t t3;
  mpz_init(t3);
  pbc_demo_pairing_init(pairing, argc, argv);
  if (!pairing_is_symmetric(pairing)) pbc_die("pairing must be symmetric");

  element_init_G1(P, pairing);
  element_init_G1(Ppub, pairing);
  element_init_G1(Qid, pairing);
  element_init_G1(Did, pairing);
  element_init_G1(R, pairing);
  element_init_G1(S, pairing);
  element_init_G1(t2, pairing);
  element_init_G1(t4, pairing);
  element_init_G1(t5, pairing);
  element_init_G1(t7, pairing);

  element_init_Zr(s, pairing);
  element_init_Zr(k, pairing);
  element_init_Zr(t1, pairing);

  element_init_GT(t6, pairing);
  element_init_GT(t8, pairing);
  element_init_GT(t9, pairing);
  element_init_GT(t10, pairing);
  element_init_GT(t11, pairing);

  time1 = pbc_get_time();
  printf("Paterson ID-based signature.\n");
  printf("KEYGEN\n");
  element_random(P);
  element_random(s);
  element_mul_zn(Ppub, P, s);
  element_printf("P = %B\n", P);
  element_printf("Ppub = %B\n", Ppub);
  element_from_hash(Qid, "ID", 2);
  element_printf("Qid = %B\n", Qid);
  element_mul_zn(Did, Qid, s);

  printf("SIGN\n");
  element_random(k);
  element_mul_zn(R, P, k);
  element_from_hash(t1, "Message", 7);
  element_mul_zn(t2, P, t1);
  //H3(R)=mpz(R);
//  int n = element_length_in_bytes(R);
//  unsigned char *data=malloc(n);
//  element_to_bytes(data, R);
//  printf("data = %s\n", data);
  element_to_mpz(t3, R);
  element_mul_mpz(t4, Did, t3);
  element_add(t5, t4, t2);
  element_invert(k, k);
  element_mul_zn(S, t5, k);
  printf("Signature of message \"Message\" is:\n");
  element_printf("R = %B\n", R);
  element_printf("S = %B\n", S);

  printf("VERIFY\n");
  element_from_hash(t1, "Message", 7);
  element_mul_zn(t7, P, t1);
  element_pairing(t6, P, t7);
  element_pairing(t8, Ppub, Qid);
  element_to_mpz(t3, R);
  element_pow_mpz(t9, t8, t3);
  element_printf("t8 = %B\n", t8);
  element_printf("t9 = %B\n", t9);
  element_mul(t10, t6, t9);
  element_printf("t10 = %B\n", t10);
  element_pairing(t11, R, S);
  element_printf("[e(P, P)^H2(M)][e(Ppub, Qid)^H3(R)] = %B\n", t10);
  element_printf("e(R, S) = %B\n", t11);
  if (!element_cmp(t10, t11)) {
    printf("Signature is valid!\n");
  } else {
    printf("Signature is invalid!\n");
  }
  time2 = pbc_get_time();
  printf("All time = %fs\n", time2 - time1);

  element_clear(P);
  element_clear(Ppub);
  element_clear(Qid);
  element_clear(Did);
  element_clear(R);
  element_clear(t2);
  element_clear(t4);
  element_clear(t5);
  element_clear(s);
  element_clear(k);
  element_clear(t1);
  element_clear(t6);
  element_clear(t7);
  element_clear(t8);
  element_clear(t9);
  element_clear(t10);
  element_clear(t11);
  pairing_clear(pairing);

  return 0;
}
