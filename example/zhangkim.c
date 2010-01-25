/*
  Zhang and Kim ID-based Blind Signature scheme.
  Based on papers "F. Zang, K. Kim. ID-based Blind Signature and Ring Signature from Pairings. Advances in Cryptology - Asiacrypt 2002, LNCS Vol. 2510, Springer-Verlag, 2002".
  Contributed by Dmitry Kosolapov.
*/

#include <pbc.h>
#include <pbc_test.h>

int main(int argc, char **argv) {
  pairing_t pairing;
  double time1, time2;
  element_t P, Ppub, s, R, Qid, Sid, a, b, r, c, S, negc, t1, t2, t3, t5,
    t6, t7, t8, t9, t10, t11, t12, t14;
  mpz_t t4, t13;
  mpz_init(t4);
  mpz_init(t13);
  pbc_demo_pairing_init(pairing, argc, argv);
  if (!pairing_is_symmetric(pairing)) pbc_die("pairing must be symmetric");
  element_init_G1(P, pairing);
  element_init_G1(Ppub, pairing);
  element_init_G1(Qid, pairing);
  element_init_G1(Sid, pairing);
  element_init_G1(R, pairing);
  element_init_G1(S, pairing);
  element_init_G1(t1, pairing);
  element_init_G1(t2, pairing);
  element_init_G1(t7, pairing);
  element_init_G1(t8, pairing);
  element_init_G1(t9, pairing);

  element_init_Zr(r, pairing);
  element_init_Zr(s, pairing);
  element_init_Zr(c, pairing);
  element_init_Zr(a, pairing);
  element_init_Zr(b, pairing);
  element_init_Zr(negc, pairing);
  element_init_Zr(t5, pairing);
  element_init_Zr(t6, pairing);
  element_init_Zr(t14, pairing);

  element_init_GT(t3, pairing);
  element_init_GT(t10, pairing);
  element_init_GT(t11, pairing);
  element_init_GT(t12, pairing);

  time1 = pbc_get_time();
  printf("Zhang and Kim ID-based Blind Signature scheme\n");
  printf("SETUP\n");
  element_random(P);
  element_random(s);
  element_mul_zn(Ppub, P, s);
  element_printf("P = %B\n", P);
  element_printf("Ppub = %B\n", Ppub);

  printf("EXTRACT\n");
  element_from_hash(Qid, "ID", 2);
  element_mul_zn(Sid, Qid, s);
  element_printf("Public key Qid = %B\n", Qid);
  element_printf("Private key Sid = %B\n", Sid);

  printf("BLIND SIGNATURE ISSUING PROTOCOL\n");
  element_random(r);
  element_mul_zn(R, P, r);
  printf("Signer sends R = rP to user\n");
  element_printf("R = %B\n", R);
  printf("Blinding\n");
  element_random(a);
  element_random(b);
  element_mul_zn(t1, P, a);
  element_add(t1, R, t1);
  element_mul_zn(t2, Qid, b);
  element_add(t2, t2, t1);
  element_pairing(t3, t2, Ppub);
  element_to_mpz(t4, t3);
  element_from_hash(t5, "Message", 7);
  element_mul_mpz(t6, t5, t4);
  element_add(c, t6, b);
  printf("User sends c to signer\n");
  element_printf("c = %B\n", c);
  printf("Signing\n");
  element_mul_zn(t7, Ppub, r);
  element_mul_zn(t8, Sid, c);
  element_add(S, t8, t7);
  printf("Signer sends S\n");
  element_printf("S = %B\n", S);
  printf("Unblinding\n");
  element_mul_zn(t9, Ppub, a);
  element_add(S, S, t9);
  element_sub(c, c, b);
  printf("Blind Signature of message \"Message\" is:\n");
  element_printf("S1 = %B\n", S);
  element_printf("c1 = %B\n", c);

  printf("VERIFICATION\n");
  element_pairing(t10, Qid, Ppub);
  element_neg(negc, c);
  element_pow_zn(t10, t10, negc);
  element_pairing(t11, S, P);
  element_mul(t12, t11, t10);
  element_to_mpz(t13, t12);
  element_from_hash(t5, "Message", 7);
  element_mul_mpz(t14, t5, t13);
  element_printf("c1 = %B\n", c);
  element_printf("H(m, [e(S1, P)][e(Qid, Ppub)^(-c1)]) = %B\n", t14);

  if (!element_cmp(t14, c)) printf("Signature is valid\n");
  else printf("Signature is invalid\n");
  time2 = pbc_get_time();
  printf("All time = %fs\n", time2 - time1);

  element_clear(P);
  element_clear(Ppub);
  element_clear(Qid);
  element_clear(Sid);
  element_clear(R);
  element_clear(S);
  element_clear(r);
  element_clear(s);
  element_clear(c);
  element_clear(a);
  element_clear(b);
  element_clear(negc);
  element_clear(t1);
  element_clear(t2);
  element_clear(t3);
  element_clear(t5);
  element_clear(t6);
  element_clear(t14);
  element_clear(t7);
  element_clear(t8);
  element_clear(t9);
  element_clear(t10);
  element_clear(t11);
  element_clear(t12);
  pairing_clear(pairing);

  return 0;
}
