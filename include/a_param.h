#ifndef A_PARAM_H
#define A_PARAM_H

#include "fieldquadratic.h"
#include "pairing.h"
#include "fops.h"

struct a_param_s {
    int exp2;
    int exp1;
    int sign1;
    int sign0;
    mpz_t r; // r = 2^exp2 + sign1 * 2^exp1 + sign0 * 1
    mpz_t q; // we work in E(F_q) (and E(F_q^2))
    mpz_t h; // r * h = q + 1
};
typedef struct a_param_s a_param_t[1];
typedef struct a_param_s *a_param_ptr;

/*@manual aparam
Initialize ''p''. This must be called before ''p'' can be used.
*/
void a_param_init(a_param_t p);

/*@manual aparam
Clear ''p''. This should be called after ''p'' is no longer needed.
*/
void a_param_clear(a_param_t p);

/*@manual aparam
Generate type A pairing parameters and store them in ''p'',
where the group order r is ''rbits'' long, and the order of the base field q
is ''qbits'' long. To be secure, generic discrete log algorithms must
be infeasible in groups of order r, and finite field discrete log algorithms
must be infeasible in finite fields of order q^2.
Typical values: ''rbits'' = 160, ''qbits'' = 512.
*/
void a_param_gen(a_param_t p, int rbits, int qbits);

/*@manual aparam
Write the parameters in ''p'' in a text format onto ''stream''.
*/
void a_param_out_str(FILE *stream, a_param_ptr p);

/*@manual aparam_internal
TODO
*/
void a_param_inp_generic (a_param_ptr p, fetch_ops_t fops, void *ctx);

/*@manual aparam_internal
Initializes ''pairing'' with type A parameters in ''p''.
*/
void pairing_init_a_param(pairing_t pairing, a_param_t p);

#endif //A_PARAM_H
