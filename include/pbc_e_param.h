// Type E pairings.

// Requires:
// * fops.h
// * param.h
#ifndef __PBC_E_PARAM_H__
#define __PBC_E_PARAM_H__

void pbc_param_init_e(pbc_param_ptr p, fetch_ops_t fops, void *ctx);

/*@manual eparam
Generate type E pairing parameters and store them in 'p',
where the group order r is 'rbits' long, and the order of the base field q
is 'qbits' long. To be secure, generic discrete log algorithms must
be infeasible in groups of order r, and finite field discrete log algorithms
must be infeasible in finite fields of order q.

Typical values: 'rbits' = 160, 'qbits' = 1024.
*/
void pbc_param_init_e_gen(pbc_param_t p, int rbits, int qbits);

#endif //__PBC_E_PARAM_H__
