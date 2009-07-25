// Type F pairings.

// Requires:
// * param.h
#ifndef __PBC_F_PARAM_H__
#define __PBC_F_PARAM_H__

void pbc_param_init_f(pbc_param_ptr par, const char *(*tab)(const char *));

/*@manual fparam
Generate type F pairing parameters and store them in 'p'.
Both the group order r and the order of the base field q will be roughly
'bits'-bit numbers.
To be secure, generic discrete log algorithms must
be infeasible in groups of order r, and finite field discrete log algorithms
must be infeasible in finite fields of order q^12.
Typical value: 'bits' = 160.
*/
void pbc_param_init_f_gen(pbc_param_t p, int bits);

#endif //__PBC_F_PARAM_H__
