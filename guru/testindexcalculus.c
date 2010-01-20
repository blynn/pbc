#include <stdio.h>
#include <gmp.h>
#include "pbc.h"

int main(int argc, char **argv)
{
    mpz_t x;
    mpz_t g, h, q;
    mpz_init(x);
    mpz_init(g);
    mpz_init(h);
    mpz_init(q);
    int bits = 40;

    if (argc == 2) {
        bits = atoi(argv[1]);
    }
    mpz_setbit(q, bits);
    pbc_mpz_random(q, q);
    mpz_nextprime(q, q);
    pbc_mpz_random(g, q);
    pbc_mpz_random(h, q);
    mpz_powm(h, g, h, q);

    element_dlog_index_calculus(x, g, h, q);
    element_printf("%Zd^%Zd %% %Zd = %Zd\n", g, x, q, h);

    return 0;
}
