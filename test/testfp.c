#include "pbc.h"
#include "fp.h"
#include "get_time.h"

void timefield(field_t fp)
{
    int i, n;
    double t0, t1;

    element_t x, y, z;
    element_init(x, fp);
    element_init(y, fp);
    element_init(z, fp);

    element_random(x);
    element_random(y);

    n = 20000;
    t0 = get_time();
    for (i=0; i<n; i++) {
	element_mul(z, x, y);
	element_mul(x, y, z);
	element_mul(y, z, x);
    }
    t1 = get_time();
    printf("mul %fs\n", t1 - t0);

    n = 20000;
    t0 = get_time();
    for (i=0; i<n; i++) {
	element_square(x, x);
    }
    t1 = get_time();
    printf("square %fs\n", t1 - t0);

    /*
    n = 10000;
    t0 = get_time();
    for (i=0; i<n; i++) {
	element_invert(z, x);
	element_invert(z, y);
    }
    t1 = get_time();
    printf("invert %fs\n", t1 - t0);
    */

    n = 40000;
    t0 = get_time();
    for (i=0; i<n; i++) {
	element_set(z, x);
	element_set(z, y);
    }
    t1 = get_time();
    printf("set %fs\n", t1 - t0);

    element_clear(x);
    element_clear(y);
    element_clear(z);
}

int main(void)
{
    field_t f1, f2;
    mpz_t prime;

    mpz_init(prime);
    mpz_setbit(prime, 1023);
    mpz_setbit(prime, 70);
    mpz_nextprime(prime, prime);
    field_init_fast_fp(f1, prime);
    field_init_naive_fp(f2, prime);

    printf("fastfp.c\n");
    timefield(f1);
    printf("naivefp.c\n");
    timefield(f2);
    //mem_report();
    return 0;
}
