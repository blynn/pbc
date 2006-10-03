#include <pbc.h>
#include "fp.h"
#include "random.h"
#include "get_time.h"

int main(void)
{
    mpz_t p, q, N, z;
    field_t f;
    element_t a;
    double t0, t1, ttotal = 0;
    int i, n;

    mpz_init(p);
    mpz_init(q);
    mpz_init(N);
    mpz_init(z);
    mpz_setbit(N, 513);
    mpz_setbit(N, 513);
    pbc_mpz_random(p, N);
    pbc_mpz_random(q, N);
    mpz_nextprime(p, p);
    mpz_nextprime(q, q);
    mpz_mul(N, p, q);

    field_init_fp(f, N);
    element_init(a, f);
    //field_print_info(stdout, f);
    n = 10;
    for (i=0; i<n; i++) {
	pbc_mpz_random(z, N);
	element_random(a);
	t0 = get_time();
	element_pow_mpz(a, a, z);
	t1 = get_time();
	ttotal += t1 - t0;
    }
    printf("average RSA exp time = %lf\n", ttotal / n);
    return 0;
}
