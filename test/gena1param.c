#include "pbc.h"
#include "a1_param.h"

int main(void)
{
    gmp_randstate_t rstate;
    a1_param_t param;
    mpz_t p, q, N;

    mpz_init(p);
    mpz_init(q);
    mpz_init(N);

    gmp_randinit_default(rstate);
    mpz_urandomb(p, rstate, 512);
    mpz_nextprime(p, p);
    mpz_urandomb(q, rstate, 512);
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
