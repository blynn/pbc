#include "pbc.h"

int main(void)
{
    mpz_t p, q, N;

    mpz_init(p);
    mpz_init(q);
    mpz_init(N);

    // In a real application, p and q must be stored somewhere safe.
    pbc_mpz_randomb(p, 512);
    pbc_mpz_randomb(q, 512);

    mpz_nextprime(p, p);
    mpz_nextprime(q, q);
    mpz_mul(N, p, q);

    pbc_param_t param;
    pbc_param_init_a1_gen(param, N);
    pbc_param_out_str(stdout, param);
    pbc_param_clear(param);
    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(N);
    return 0;
}
