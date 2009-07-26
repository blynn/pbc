// Test F_p.

#include "pbc.h"
#include "pbc_fp.h"
#include "pbc_test.h"

int main(void) {
  field_t fp;
  mpz_t prime;
  mpz_t m, n;

  mpz_init(prime);
  mpz_init(n);
  mpz_init(m);
  mpz_set_ui(prime, 100000);
  mpz_setbit(prime, 33);
  mpz_nextprime(prime, prime);

  field_init_fp(fp, prime);

  element_t x, y, z;
  element_init(x, fp);
  element_init(y, fp);
  element_init(z, fp);

  long a = 123, b = 456;

  // Conversion to and from signed long.
  EXPECT(0 == element_to_si(z));
  element_set1(z);
  EXPECT(1 == element_to_si(z));
  element_set0(z);
  EXPECT(0 == element_to_si(z));
  element_set_si(x, a);
  EXPECT(a == element_to_si(x));
  element_set_si(y, b);
  EXPECT(b == element_to_si(y));
  // Assignment, comparison.
  EXPECT(!element_cmp(x, x));
  EXPECT(element_cmp(x, y));
  EXPECT(element_cmp(z, x));
  element_set(z, x);
  EXPECT(!element_cmp(z, x));
  // Arithmetic operations.
  element_add(z, x, y);
  EXPECT(a + b == element_to_si(z));
  element_mul(z, x, y);
  EXPECT(a * b == element_to_si(z));
  element_sub(z, y, x);
  EXPECT(b - a == element_to_si(z));
  element_set_mpz(z, prime);
  EXPECT(!element_to_si(z));
  element_sub(z, z, x);
  element_to_mpz(n, z);
  mpz_add_ui(n, n, a);
  EXPECT(!mpz_cmp(n, prime));
  element_invert(z, x);
  element_to_mpz(m, z);
  mpz_set_ui(n, a);
  mpz_invert(n, n, prime);
  EXPECT(!mpz_cmp(m, n));
  element_invert(z, z);
  EXPECT(!element_cmp(x, z));
  element_div(z, y, x);
  element_to_mpz(m, z);
  mpz_mul_ui(n, n, b);
  mpz_mod(n, n, prime);
  EXPECT(!mpz_cmp(m, n));
  // Exponentiation.
  element_pow_zn(z, x, y);
  element_to_mpz(m, z);
  mpz_set_si(n, a);
  mpz_powm_ui(n, n, b, prime);
  EXPECT(!mpz_cmp(m, n));
  // Preprocessed exponentiation.
  element_pp_t p;
  element_pp_init(p, x);
  element_pp_pow_zn(z, y, p);
  element_pp_clear(p);
  element_to_mpz(m, z);
  EXPECT(!mpz_cmp(m, n));

  element_from_hash(z, NULL, 0);
  element_from_hash(x, NULL, 0);
  EXPECT(!element_cmp(z, x));

  element_clear(x);
  element_clear(y);
  element_clear(z);
  field_clear(fp);
  mpz_clear(prime);
  mpz_clear(m);
  mpz_clear(n);
  return pbc_err_count;
}
