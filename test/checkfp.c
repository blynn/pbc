//compares two implementations of Fp against each other
#include <string.h>
#include "pbc.h"
#include "fp.h"
#include "fieldquadratic.h"

static mpz_t prime;

static void check_match(element_t e1, element_t e2, char *s)
{
    unsigned char *buf1, *buf2;
    int len;
    void bug(void)
    {
	printf("BUG: %s mismatch\n", s);
	element_printf("e1: %B\n", e1);
	element_printf("e2: %B\n", e2);
	exit(1);
    }

    len = element_length_in_bytes(e1);
    if (len != element_length_in_bytes(e2)) {
	bug();
    }

    buf1 = malloc(len);
    buf2 = malloc(len);
    element_to_bytes(buf1, e1);
    element_to_bytes(buf2, e2);

    if (memcmp(buf1, buf2, len)) {
	bug();
    }

    free(buf1);
    free(buf2);
}

static void run_check(field_ptr f1, field_ptr f2)
{
    mpz_t t1, t2;
    element_t x1, y1, z1;
    element_t x2, y2, z2;
    unsigned char *buf;
    int len;

    mpz_init(t1);
    mpz_init(t2);
    element_init(x1, f1);
    element_init(y1, f1);
    element_init(z1, f1);
    element_init(x2, f2);
    element_init(y2, f2);
    element_init(z2, f2);

    check_match(z1, z2, "init");
    element_random(x1);
    len = element_length_in_bytes(x1);
    buf = malloc(len);
    element_to_bytes(buf, x1);
    element_from_bytes(x2, buf);
    free(buf);
    check_match(x1, x2, "conversion");
    element_random(y1);
    len = element_length_in_bytes(y1);
    buf = malloc(len);
    element_to_bytes(buf, y1);
    element_from_bytes(y2, buf);
    free(buf);
    check_match(y1, y2, "conversion");
    element_add(z1, x1, y1);
    element_add(z2, x2, y2);
    check_match(z1, z2, "add");
    element_add(z1, z1, z1);
    element_add(z2, z2, z2);
    check_match(z1, z2, "add (in place)");
    element_sub(z1, x1, y1);
    element_sub(z2, x2, y2);
    check_match(z1, z2, "sub");
    element_mul(z1, x1, y1);
    element_mul(z2, x2, y2);
    check_match(z1, z2, "mul");
    element_mul(z1, z1, z1);
    element_mul(z2, z2, z2);
    check_match(z1, z2, "mul (in place)");
    element_neg(y1, x1);
    element_neg(y2, x2);
    check_match(y1, y2, "neg");
    element_add(z1, x1, y1);
    element_add(z2, x2, y2);
    check_match(z1, z2, "add (to zero)");
    element_invert(y1, x1);
    element_invert(y2, x2);
    check_match(y1, y2, "invert");
    element_mul(z1, x1, y1);
    element_mul(z2, x2, y2);
    check_match(z1, z2, "mul (to one)");
    element_square(z1, x1);
    element_square(z2, x2);
    check_match(z1, z2, "square");
    element_double(z1, x1);
    element_double(z2, x2);
    check_match(z1, z2, "double");
    while (!element_is_sqr(x1)) {
	element_random(x1);
    }
    len = element_length_in_bytes(x1);
    buf = malloc(len);
    element_to_bytes(buf, x1);
    element_from_bytes(x2, buf);
    free(buf);
    check_match(x1, x2, "conversion");
    element_sqrt(z1, x1);
    element_sqrt(z2, x2);
    check_match(z1, z2, "sqrt");
    element_to_mpz(t1, y1);
    element_to_mpz(t2, y2);
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
    
int main(void)
{
    field_t f1, f2;
    field_t f1i, f2i;
    field_t f1x, f2x;
    field_t f1p, f2p;
    int i, n;
    element_ptr n1;
    element_t n2;
    element_t irred1, irred2;
    mpz_t z;

    mpz_init(z);
    mpz_init(prime);
    mpz_set_ui(prime, 82);
    mpz_setbit(prime, 32);
    mpz_nextprime(prime, prime);

    field_init_naive_fp(f1, prime);
    field_init_fast_fp(f2, prime);

    element_init(n2, f2);

    n1 = field_get_nqr(f1);
    element_to_mpz(z, n1);
    element_set_mpz(n2, z);
    field_set_nqr(f2, n2);

    field_init_fi(f1i, f1);
    field_init_fi(f2i, f2);

    n = 10;

    element_printf("prime = %Z\n", prime);
    printf("checking base fields\n");
    for (i=0; i<n; i++) run_check(f1, f2);

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
	buf = malloc(len);
	element_to_bytes(buf, irred1);
	element_from_bytes(irred2, buf);
	free(buf);
    }
    field_init_polymod(f2p, irred2);

    run_check(f1p, f2p);

    return 0;
}
