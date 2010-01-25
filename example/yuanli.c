/* Contributed by Dmitry Kosolapov
 *
 * I haven't tested this much, and I'm personally not familiar with
 * this particular cryptosystem. -Ben Lynn
 */
/* Here we represent the original Yuan-Li ID-Based Authenticated Key Agreement Protocol, 2005.
 * This protocol has 2 stages: Setup and Extract. We represent them inside one code block with demo and time outputs.
 */

/*Yuan-Li protocol description according to:
Quan Yuan and Songping Li, A New Efficient ID-Based Authenticated Key Agreement Protocol, Cryptology ePrint Archive, Report 2005/309

SETUP:
KGS chooses G1, G2, e: G1*G1 -> G2, P, H: {0, 1}* -> G1, s, H - some function for key calculation.
KGS calculates Ppub = s*P, publishes {G1, G2, e, P, Ppub, H1, H} and saves s as master key.

EXTRACT:

For the user with ID public key can be calculated with Qid = H1(ID). KGS generates bound public key Sid = s*Qid.
1. A chooses random a from Z_p*, calculates Ta = a*P.
   A -> B: Ta
2. B chooses random b from Z_p*, calculates Tb = b*P.
   B -> A: Tb
3. A calculates h = a*Tb = a*b*P and shared secret key Kab = e(a*Ppub + Sa, Tb + Qb)
4. B calculates h = b*Ta = a*b*P and shared secret key Kba = e(Ta + Qa, b*Ppub + Sb)
Session key is K = H(A, B, h, Kab).
H was not defined in the original article.
I've defined it as H(A, B, h, Kab)=e(h,H1(A)+H1(B))+Kab.
*/

#include <pbc.h>
#include <pbc_test.h>

int main(int argc, char **argv) {
  pairing_t pairing;
  double t0, t1;
  element_t s, a, b, P, Ppub, Qa, Qb, Sa, Sb, Ta, Tb, Kab, Kba, K, temp1,
    temp2, temp3, temp4, temp5, h;

  pbc_demo_pairing_init(pairing, argc, argv);
  if (!pairing_is_symmetric(pairing)) pbc_die("pairing must be symmetric");

  element_init_Zr(s, pairing);
  element_init_Zr(a, pairing);
  element_init_Zr(b, pairing);

  element_init_G1(P, pairing);
  element_init_G1(Ppub, pairing);
  element_init_G1(Qa, pairing);
  element_init_G1(Qb, pairing);
  element_init_G1(Sa, pairing);
  element_init_G1(Sb, pairing);
  element_init_G1(Ta, pairing);
  element_init_G1(Tb, pairing);
  element_init_G1(temp1, pairing);
  element_init_G1(temp2, pairing);
  element_init_G1(temp3, pairing);
  element_init_G1(h, pairing);

  element_init_GT(Kab, pairing);
  element_init_GT(Kba, pairing);
  element_init_GT(K, pairing);
  element_init_GT(temp4, pairing);
  element_init_GT(temp5, pairing);

  printf("Yuan-Li key agreement protocol\n");

  t0 = pbc_get_time();

//Setup, system parameters generation
  printf("SETUP STAGE\n");
  element_random(P);
  element_printf("P = %B\n", P);
  element_random(s);
  element_mul_zn(Ppub, P, s);
  element_printf("Ppub = %B\n", Ppub);

//Extract, key calculation
  printf("EXTRACT STAGE\n");
  element_from_hash(Qa, "A", 1);
  element_from_hash(Qb, "B", 1);
  element_mul_zn(Sa, Qa, s);
  element_mul_zn(Sb, Qb, s);
  element_printf("Sa = %B\n", Sa);
  element_printf("Sb = %B\n", Sb);

  printf("-----1-----\n");

  element_random(a);
  element_mul_zn(Ta, P, a);
  element_printf("A sends B Ta = %B\n", Ta);

  printf("-----2-----\n");

  element_random(b);
  element_mul_zn(Tb, P, b);
  element_printf("B sends A Tb = %B\n", Tb);

  printf("-----3-----\n");

  printf("A calculates h and Kab\n");
  element_mul_zn(h, Tb, a);
  element_printf("h = %B\n", h);
  element_mul_zn(temp1, Ppub, a);
  element_add(temp1, temp1, Sa);
  element_add(temp2, Tb, Qb);
  element_pairing(Kab, temp1, temp2);
  element_printf("Kab = %B\n", Kab);

  printf("-----4-----\n");

  printf("B calculates h and Kba\n");
  element_mul_zn(h, Ta, b);
  element_printf("h = %B\n", h);
  element_add(temp1, Ta, Qa);
  element_mul_zn(temp2, Ppub, b);
  element_add(temp2, temp2, Sb);
  element_pairing(Kba, temp1, temp2);
  element_printf("Kba = %B\n", Kba);

  printf("-----FINAL-----\n");

  element_add(temp3, Qa, Qb);
  element_pairing(temp4, h, temp3);

  element_add(K, temp4, Kab);
  element_printf("A has the key K = %B\n", K);
  element_set(temp5, K);

  element_add(K, temp4, Kba);
  element_printf("B has the key K = %B\n", K);

  if (!element_cmp(temp5, K))
    printf("The keys are the same. Start session...\n");
  else
    printf("The keys aren't the same. Try again, please.\n");

  element_clear(K);
  element_clear(Kab);
  element_clear(Kba);
  element_clear(h);
  element_clear(temp1);
  element_clear(temp2);
  element_clear(temp3);
  element_clear(temp4);
  element_clear(temp5);
  element_clear(s);
  element_clear(a);
  element_clear(b);
  element_clear(P);
  element_clear(Ppub);
  element_clear(Qa);
  element_clear(Qb);
  element_clear(Sa);
  element_clear(Sb);
  element_clear(Ta);
  element_clear(Tb);

  t1 = pbc_get_time();

  printf("All time = %fs\n", t1 - t0);
  printf("Have a good day!\n");

  return 0;
}
