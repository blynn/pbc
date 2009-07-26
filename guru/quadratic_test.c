// Test quadratic field extensions.

#include "pbc.h"
#include "pbc_fp.h"
#include "pbc_fieldquadratic.h"
#include "pbc_test.h"

int main(void) {
  field_t fp, fp2;
  mpz_t prime;
  element_t a, b, c;

  mpz_init(prime);
  // Prime is 3 mod 4 so that -1 is a quadratic nonresidue.
  // For smaller tests, try the prime 83.
  mpz_setbit(prime, 256);
  do {
    mpz_nextprime(prime, prime);
  } while (mpz_fdiv_ui(prime, 4) != 3);

  field_init_fp(fp, prime);
  field_init_fi(fp2, fp);
  element_init(a, fp2);
  element_init(b, fp2);
  element_init(c, fp2);

  element_printf("field: %Z^2\n", prime);

  element_random(a);
  element_random(b);
  element_printf("a = %B, b = %B\n", a, b);

  element_add(c, a, b);
  element_printf("a + b = %B\n", c);

  element_mul(c, a, b);
  element_printf("a * b = %B\n", c);

  for (;;) {
    element_random(a);
    element_printf("new a = %B\n", a);

    if (element_is_sqr(a)) break;
    printf(" is not a square\n");
  }
  element_sqrt(c, a);
  element_printf("sqrt(a) = %B\n", c);
  element_mul(c, c, c);
  element_printf("sqrt(a) * sqrt(a) = %B\n", c);
  element_invert(c, a);
  element_printf("1/a = %B\n", c);
  element_mul(c, c, a);
  element_printf("1/a * a = %B\n", c);

  element_clear(a);
  element_clear(b);
  element_clear(c);
  field_clear(fp);
  field_clear(fp2);
  mpz_clear(prime);
  return 0;
}
