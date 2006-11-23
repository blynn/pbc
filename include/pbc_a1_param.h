// requires
// * stdio.h
// * gmp.h
// * pairing.h
// * fops.h
#ifndef __PBC_A1_PARAM_H__
#define __PBC_A1_PARAM_H__

struct a1_param_s {
    mpz_t p;
    mpz_t n;
    unsigned int l;
};
typedef struct a1_param_s a1_param_t[1];
typedef struct a1_param_s *a1_param_ptr;

/*@manual a1param
Initialize ''p''. This must be called before ''p'' can be used.
*/
void a1_param_init(a1_param_t param);

/*@manual a1param
Clear ''p''. This should be called after ''p'' is no longer needed.
*/
void a1_param_clear(a1_param_t param);

/*@manual a1param
Write the parameters in ''p'' in a text format onto ''stream''.
*/
void a1_param_out_str(FILE *stream, a1_param_ptr p);

/*@manual a1param
Generate type A1 pairing parameters and store them in ''p''.
The group order will be ''n''. The order of the base field is a few bits longer.
To be secure, generic discrete log algorithms must
be infeasible in groups of order ''n'',
and finite field discrete log algorithms
must be infeasible in finite fields of order roughly ''n''^2.
Furthermore, ''n'' should be hard to factorize.
Typical values: ''n'' is a product of two primes, each at least 512 bits long.
*/
void a1_param_gen(a1_param_t param, mpz_t n);

void a1_param_inp_generic(a1_param_ptr p, fetch_ops_t fops, void *ctx);
void pairing_init_a1_param(pairing_t pairing, a1_param_t param);

#endif //__PBC_A1_PARAM_H__
