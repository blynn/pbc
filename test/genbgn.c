#include "pbc.h"
#include "bgn_param.h"

int main(void)
{
    gmp_randstate_t rstate;
    bgn_param_t param;
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

    bgn_param_init(param);
    bgn_param_gen(param, N);

    bgn_param_out_str(stdout, param);

    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(N);
    bgn_param_clear(param);
    return 0;
}
