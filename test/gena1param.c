#include "pbc.h"

int main(void)
{
    a1_param_t param;
    mpz_t p, q, N;

    mpz_init(p);
    mpz_init(q);
    mpz_init(N);

    //in a real application, p and q must be stored somewhere safe
    pbc_mpz_randomb(p, 512);
    pbc_mpz_randomb(q, 512);

    mpz_nextprime(p, p);
    mpz_nextprime(q, q);
    mpz_mul(N, p, q);

    a1_param_init(param);
    a1_param_gen(param, N);

    a1_param_out_str(stdout, param);

    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(N);
    a1_param_clear(param);
    return 0;
}
