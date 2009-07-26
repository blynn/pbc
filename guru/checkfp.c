// Compares two implementations of Fp.

#include <string.h>
#include "pbc.h"
#include "pbc_fp.h"
#include "pbc_fieldquadratic.h"

static mpz_t prime;

enum { VERBOSE = 0 };

static void check_p(int value, char *s) {
  if (!value) {
    printf("BUG: %s predicate wrong\n", s);
    exit(1);
  }

  if (VERBOSE) {
    printf("checking %s\n", s);
  }
}

static void check_match_int(int i1, int i2, char *s) {
  void bug(void)
  {
    printf("BUG: %s mismatch\n", s);
    element_printf("i1: %d\n", i1);
    element_printf("i2: %d\n", i2);
    exit(1);
  }

  if (VERBOSE) {
    printf("checking %s\n", s);
    element_printf("i1: %d\n", i1);
    element_printf("i2: %d\n", i2);
  }

  if (i1 != i2) bug();
}

static void check_match(element_t e1, element_t e2, char *s) {
  unsigned char *buf1, *buf2;
  int len;
  void bug(void)
  {
    printf("BUG: %s mismatch\n", s);
    element_printf("e1: %B\n", e1);
    element_printf("e2: %B\n", e2);
    exit(1);
  }

  if (VERBOSE) {
    printf("checking %s\n", s);
    element_printf("e1: %B\n", e1);
    element_printf("e2: %B\n", e2);
  }
  len = element_length_in_bytes(e1);
  if (len != element_length_in_bytes(e2)) {
    bug();
  }

  buf1 = pbc_malloc(len);
  buf2 = pbc_malloc(len);
  element_to_bytes(buf1, e1);
  element_to_bytes(buf2, e2);

  if (memcmp(buf1, buf2, len)) {
    bug();
  }

  pbc_free(buf1);
  pbc_free(buf2);
}

static void run_check(field_ptr f1, field_ptr f2) {
  mpz_t t1, t2;
  element_t x1, y1, z1;
  element_t x2, y2, z2;
  char s2[80];

  void convertset(element_t out, element_t in)
  {
    unsigned char *buf;
    int len;

    len = element_length_in_bytes(in);
    buf = pbc_malloc(len);
    element_to_bytes(buf, in);
    element_from_bytes(out, buf);
    pbc_free(buf);
    check_match(in, out, "conversion");
  }

  void randxy(void)
  {

    element_random(x1);
    element_random(y1);
    convertset(x2, x1);
    convertset(y2, y1);
  }

  void check_onearg(void (*fn)(element_ptr), char *s)
  {
    fn(x1);
    fn(x2);
    check_match(x1, x2, s);
  }

  void check_twoarg(void (*fn)(element_ptr, element_ptr), char *s)
  {
    randxy();
    fn(z1, x1);
    fn(z2, x2);
    check_match(z1, z2, s);

    strncpy(s2, s, 32);
    strcat(s2, " (in place)");
    fn(y1, y1);
    fn(y2, y2);
    check_match(y1, y2, s2);
  }

  void check_threearg(void (*fn)(element_ptr, element_ptr, element_ptr), char *s)
  {
    randxy();
    fn(z1, x1, y1);
    fn(z2, x2, y2);
    check_match(z1, z2, s);

    strncpy(s2, s, 32);
    strcat(s2, " (first arg in place)");
    element_set(z1, x1);
    element_set(z2, x2);
    fn(z1, z1, y1);
    fn(z2, z2, y2);
    check_match(z1, z2, s2);

    strncpy(s2, s, 32);
    strcat(s2, " (second arg in place)");
    element_set(z1, y1);
    element_set(z2, y2);
    fn(z1, x1, z1);
    fn(z2, x2, z2);
    check_match(z1, z2, s2);

    strncpy(s2, s, 32);
    strcat(s2, " (both args in place)");
    element_set(z1, y1);
    element_set(z2, y2);
    fn(x1, x1, x1);
    fn(x2, x2, x2);
    check_match(x1, x2, s2);
  }

  mpz_init(t1);
  mpz_init(t2);
  element_init(x1, f1);
  element_init(y1, f1);
  element_init(z1, f1);
  element_init(x2, f2);
  element_init(y2, f2);
  element_init(z2, f2);

  check_p(!element_cmp(x1, y1), "cmp0-1");
  check_p(!element_cmp(x2, y2), "cmp0-2");
  check_match(z1, z2, "init");
  check_onearg(element_set0, "set0");
  check_onearg(element_set1, "set1");
  check_twoarg(element_set, "set");
  check_match_int(element_sgn(z1), element_sgn(z2), "sgn");

  check_threearg(element_add, "add");
  check_twoarg(element_neg, "neg");
  check_threearg(element_sub, "sub");
  check_twoarg(element_double, "double");
  check_twoarg(element_halve, "halve");

  check_twoarg(element_invert, "invert");
  check_twoarg(element_square, "square");
  check_threearg(element_mul, "mul");

  randxy();
  element_neg(y1, x1);
  element_neg(y2, x2);
  element_add(z1, x1, y1);
  element_add(z2, x2, y2);
  check_match(z1, z2, "add (to zero)");
  check_p(!element_sgn(z1), "sgn");
  check_p(!element_sgn(z1), "sgn");
  check_p(element_is0(z2), "is0");
  check_p(element_is0(z2), "is0");

  randxy();
  element_invert(y1, x1);
  element_invert(y2, x2);
  element_mul(z1, x1, y1);
  element_mul(z2, x2, y2);
  check_match(z1, z2, "mul (to one)");
  check_p(element_is1(z1), "is1");
  check_p(element_is1(z2), "is1");

  randxy();
  check_p(!(!!element_cmp(x1, y1) ^ !!element_cmp(x2, y2)), "cmp");
  element_set(x1, y1);
  element_set(x2, y2);
  check_p(!element_cmp(x1, y1), "cmp");
  check_p(!element_cmp(x2, y2), "cmp");
  check_p(!element_cmp(x1, x1), "cmp (in place)");
  check_p(!element_cmp(x2, x2), "cmp (in place)");

  for (;;) {
    int flag;
    randxy();
    flag = element_is_sqr(x1);
    check_match_int(flag, element_is_sqr(x2), "is_sqr");
    if (flag) break;
  }
  convertset(x2, x1);
  element_sqrt(z1, x1);
  element_sqrt(z2, x2);
  //can't compare these because sqrt is nondeterministic
  //and there's no way easy way to preserve random state yet
  element_square(z1, z1);
  element_square(z2, z2);
  check_match(z1, z2, "sqrt");

  pbc_mpz_random(t1, f1->order);
  pbc_mpz_random(t2, f2->order);
  element_to_mpz(t1, y1);
  element_to_mpz(t2, y2);
  element_set_mpz(y1, t1);
  element_set_mpz(y2, t2);
  check_match(y1, y2, "set_mpz");
  element_mul_mpz(z1, x1, t1);
  element_mul_mpz(z2, x2, t2);
  check_match(z1, z2, "mul_mpz");
  element_pow_mpz(z1, x1, t1);
  element_pow_mpz(z2, x2, t2);
  check_match(z1, z2, "pow_mpz");
  element_mul_si(z1, x1, mpz_get_ui(t1));
  element_mul_si(z2, x2, mpz_get_ui(t2));
  check_match(z1, z2, "mul_si");
  element_set_si(z1, mpz_get_ui(t1));
  element_set_si(z2, mpz_get_ui(t2));
  check_match(z1, z2, "set_si");

  element_clear(x1);
  element_clear(y1);
  element_clear(z1);
  element_clear(x2);
  element_clear(y2);
  element_clear(z2);

  mpz_clear(t1);
  mpz_clear(t2);
}

int main(void) {
  field_t f1, f2;
  field_t f1i, f2i;
  field_t f1x, f2x;
  field_t f1p, f2p;
  int i, n;
  element_ptr n1;
  element_t n2;
  element_t irred1, irred2;
  mpz_t z;

  n = 10;

  mpz_init(z);
  mpz_init(prime);
  mpz_set_ui(prime, 1234);
  mpz_setbit(prime, 160);
  mpz_nextprime(prime, prime);

  element_printf("prime = %Zd\n", prime);

  field_init_naive_fp(f1, prime);
  field_init_fp(f2, prime);

  printf("Field 1:\n");
  field_out_info(stdout, f1);
  printf("Field 2:\n");
  field_out_info(stdout, f2);

  printf("checking base fields\n");
  for (i=0; i<n; i++) run_check(f1, f2);

  element_init(n2, f2);

  n1 = field_get_nqr(f1);
  element_to_mpz(z, n1);
  element_set_mpz(n2, z);
  field_set_nqr(f2, n2);

  field_init_fi(f1i, f1);
  field_init_fi(f2i, f2);

  printf("checking quadratic field extensions\n");
  for (i=0; i<n; i++) run_check(f1i, f2i);

  field_clear(f1i);
  field_clear(f2i);
  field_init_quadratic(f1i, f1);
  field_init_quadratic(f2i, f2);
  for (i=0; i<n; i++) run_check(f1i, f2i);

  printf("checking degree 3 extension\n");
  field_init_poly(f1x, f1);
  field_init_poly(f2x, f2);
  element_init(irred1, f1x);
  element_init(irred2, f2x);
  do {
    poly_random_monic(irred1, 3);
  } while (!poly_is_irred(irred1));

  field_init_polymod(f1p, irred1);
  {
    unsigned char *buf;
    int len;
    len = element_length_in_bytes(irred1);
    buf = pbc_malloc(len);
    element_to_bytes(buf, irred1);
    element_from_bytes(irred2, buf);
    pbc_free(buf);
  }
  field_init_polymod(f2p, irred2);

  for (i=0; i<n; i++) run_check(f1p, f2p);

  return 0;
}
