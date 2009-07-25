// requires
// * gmp.h
// * param.h
#ifndef __PBC_A1_PARAM_H__
#define __PBC_A1_PARAM_H__

void pbc_param_init_a1(pbc_param_ptr par, const char *s);

/*@manual a1param
Generate type A1 pairing parameters and store them in 'p'.  The group order
will be 'n'. The order of the base field is a few bits longer.  To be secure,
generic discrete log algorithms must be infeasible in groups of order 'n', and
finite field discrete log algorithms must be infeasible in finite fields of
order roughly 'n'^2^.  Additionally, 'n' should be hard to factorize.

For example: 'n' a product of two primes, each at least 512 bits.
*/
void pbc_param_init_a1_gen(pbc_param_t param, mpz_t n);

#endif //__PBC_A1_PARAM_H__
