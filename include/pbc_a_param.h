// Type A pairing parameters.

// Requires:
// * param.h
#ifndef __PBC_A_PARAM_H__
#define __PBC_A_PARAM_H__

int pbc_param_init_a(pbc_param_ptr par, const char *(*tab)(const char *));

/*@manual aparam
Generate type A pairing parameters and store them in 'p',
where the group order r is 'rbits' long, and the order of the base field q
is 'qbits' long. Elements take 'qbits' to represent.

To be secure, generic discrete log algorithms must be infeasible in groups of
order r, and finite field discrete log algorithms must be infeasible in finite
fields of order q^2, e.g. 'rbits' = 160, 'qbits' = 512.

The file `param/a.param` contains parameters for a type A pairing suitable for
cryptographic use.
*/
void pbc_param_init_a_gen(pbc_param_ptr par, int rbits, int qbits);

#endif //__PBC_A_PARAM_H__
