#include "pbc.h"
#include "get_time.h"

int main()
{
    field_t fp;
    mpz_t prime;
    element_t x, y, z;
    int i, n;
    double t0, t1;

    mpz_init(prime);
    mpz_setbit(prime, 1024);
    mpz_setbit(prime, 70);
    mpz_nextprime(prime, prime);
    //field_init_fast_fp(fp, prime);
    field_init_naive_fp(fp, prime);

    element_init(x, fp);
    element_init(y, fp);
    element_init(z, fp);

    element_random(x);
    element_random(y);
    element_printf("prime: %Z\n", prime);
    element_printf("x: %B\n", x);
    element_printf("y: %B\n", y);
    element_mul(z, x, y);
    element_printf("xy: %B\n", z);
    element_add(z, x, y);
    element_printf("x+y: %B\n", z);
    element_sub(z, x, y);
    element_printf("x-y: %B\n", z);
    element_invert(z, x);
    element_printf("x^-1: %B\n", z);

    n = 10000;
    t0 = get_time();
    for (i=0; i<n; i++) {
	element_mul(z, x, y);
	element_mul(x, y, z);
	element_mul(y, z, x);
    }
    t1 = get_time();
    printf("time %fs\n", t1 - t0);

    element_clear(x);
    element_clear(y);
    element_clear(z);
    //mem_report();
    return 0;
}
