#include "pbc.h"
#include "pbc_fp.h"
#include "pbc_test.h"

static void timefield(field_t fp) {
  int i, n;
  double t0, t1;

  element_t x, y, z;
  element_init(x, fp);
  element_init(y, fp);
  element_init(z, fp);

  element_random(x);
  element_random(y);

  n = 20000;
  t0 = pbc_get_time();
  for (i=0; i<n; i++) {
    element_mul(z, x, y);
    element_mul(x, y, z);
    element_mul(y, z, x);
  }
  t1 = pbc_get_time();
  printf("mul %fs\n", t1 - t0);

  n = 20000;
  t0 = pbc_get_time();
  for (i=0; i<n; i++) {
    element_square(x, x);
  }
  t1 = pbc_get_time();
  printf("square %fs\n", t1 - t0);

  n = 1000;
  t0 = pbc_get_time();
  for (i=0; i<n; i++) {
    element_invert(z, x);
    element_invert(z, y);
  }
  t1 = pbc_get_time();
  printf("invert %fs\n", t1 - t0);

  n = 40000;
  t0 = pbc_get_time();
  for (i=0; i<n; i++) {
    element_set0(z);
  }
  t1 = pbc_get_time();
  printf("set0 %fs\n", t1 - t0);

  n = 40000;
  t0 = pbc_get_time();
  for (i=0; i<n; i++) {
    element_set(z, x);
    element_set(z, y);
  }
  t1 = pbc_get_time();
  printf("set %fs\n", t1 - t0);

  n = 400;
  t0 = pbc_get_time();
  for (i=0; i<n; i++) {
    element_pow_zn(x, y, z);
  }
  t1 = pbc_get_time();
  printf("pow_zn %fs\n", t1 - t0);

  element_clear(x);
  element_clear(y);
  element_clear(z);
}

int main(int argc, char **argv) {
  field_t f1, f2;
  mpz_t prime;

  mpz_init(prime);
  if (argc > 1) {
    mpz_setbit(prime, atoi(argv[1]));
  } else {
    mpz_setbit(prime, 201);
  }
  mpz_setbit(prime, 70);
  mpz_nextprime(prime, prime);
  field_init_mont_fp(f1, prime);
  field_init_faster_fp(f2, prime);

  printf("montfp.c\n");
  timefield(f1);
  printf("fasterfp.c\n");
  timefield(f2);

  mpz_clear(prime);
  field_clear(f1);
  field_clear(f2);
  return 0;
}
