//GMP based complex floats
#include <stdio.h>
#include <gmp.h>
#include "mpc.h"

//(a+bi)(c+di) = ac - bd  + ((a+b)(c+d) - ac - bd)i
void mpc_mul(mpc_t res, mpc_t z0, mpc_t z1)
{
    mpf_t ac, bd, f0;
    mpf_init(ac);
    mpf_init(bd);
    mpf_init(f0);
    mpf_mul(ac, z0->a, z1->a);
    mpf_mul(bd, z0->b, z1->b);
    mpf_add(f0, z0->a, z0->b);
    mpf_add(res->b, z1->a, z1->b);
    mpf_mul(res->b, res->b, f0);
    mpf_sub(res->b, res->b, ac);
    mpf_sub(res->b, res->b, bd);
    mpf_sub(res->a, ac, bd);
    mpf_clear(f0);
    mpf_clear(ac);
    mpf_clear(bd);
}

void mpc_mul_2exp(mpc_t res, mpc_t z, unsigned long int e)
{
    mpf_mul_2exp(res->a, z->a, e);
    mpf_mul_2exp(res->b, z->b, e);
}

//(a+bi)^2 = (a-b)(a+b) + 2abi
void mpc_sqr(mpc_t res, mpc_t z)
{
    mpf_t f0, f1;
    mpf_init(f0);
    mpf_init(f1);
    mpf_add(f0, z->a, z->b);
    mpf_sub(f1, z->a, z->b);
    mpf_mul(f0, f0, f1);
    mpf_mul(f1, z->a, z->b);
    mpf_set(res->a, f0);
    mpf_add(res->b, f1, f1);
    mpf_clear(f0);
    mpf_clear(f1);
}

//1/(a+bi) = (1/(a^2 + b^2))(a-bi)
//naive. TODO: use one that is less prone to (over/under)flows/precision loss
void mpc_inv(mpc_t res, mpc_t z)
{
    mpf_t f0, f1;
    mpf_init(f0);
    mpf_init(f1);
    mpf_mul(f0, z->a, z->a);
    mpf_mul(f1, z->b, z->b);
    mpf_add(f0, f0, f1);
    mpf_ui_div(f0, 1, f0);
    mpf_mul(res->a, z->a, f0);
    mpf_neg(f0, f0);
    mpf_mul(res->b, z->b, f0);
    mpf_clear(f0);
    mpf_clear(f1);
}

void mpc_div(mpc_t res, mpc_t z0, mpc_t z1)
{
    mpc_t c0;
    mpc_init(c0);
    mpc_inv(c0, z1);
    mpc_mul(res, z0, c0);
    mpc_clear(c0);
}

size_t mpc_out_str(FILE *stream, int base, size_t n_digits, mpc_t op)
{
    size_t result, status;
    result = mpf_out_str(stream, base, n_digits, op->a);
    if (!result) return 0;
    if (mpf_sgn(op->b) >= 0) {
        if (EOF == fputc('+', stream)) return 0;
        result++;
    }
    status = mpf_out_str(stream, base, n_digits, op->b);
    if (!status) return 0;
    if (EOF == fputc('i', stream)) return 0;
    return result + status + 1;
}

void mpc_pow_ui(mpc_t res, mpc_t z, unsigned int n)
{
    unsigned int m;
    mpc_t z0;
    mpc_init(z0);

    //set m to biggest power of 2 less than n
    for (m = 1; m <= n; m <<= 1);
    m >>= 1;

    mpf_set_ui(z0->a, 1);
    mpf_set_ui(z0->b, 0);
    while (m) {
        mpc_mul(z0, z0, z0);
        if (m & n) {
            mpc_mul(z0, z0, z);
        }
        m >>= 1;
    }
    mpc_set(res, z0);
    mpc_clear(z0);
}

void mpc_muli(mpc_t res, mpc_t z)
{
    //i(a+bi) = -b + ai
    mpf_t f0;
    mpf_init(f0);
    mpf_neg(f0, z->b);
    mpf_set(res->b, z->a);
    mpf_set(res->a, f0);
    mpf_clear(f0);
}
