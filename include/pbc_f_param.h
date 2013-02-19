// Type F pairings.

// Requires:
// * param.h
#ifndef __PBC_F_PARAM_H__
#define __PBC_F_PARAM_H__

struct symtab_s;
int pbc_param_init_f(pbc_param_ptr par, struct symtab_s *tab);

/*@manual fparam
Generate type F pairing parameters and store them in 'p'.
Both the group order r and the order of the base field q will be roughly
'bits'-bit numbers.
To be secure, generic discrete log algorithms must
be infeasible in groups of order r, and finite field discrete log algorithms
must be infeasible in finite fields of order q^12, e.g. 'bits' = 160.

Type F should be used when the top priority is to minimize bandwidth (e.g.
short signatures). The current implementation makes them slow.

If finite field discrete log algorithms improve further, type D pairings will
have to use larger fields, but type F can still remain short, up to a point.
*/
void pbc_param_init_f_gen(pbc_param_t p, int bits);

#endif //__PBC_F_PARAM_H__
