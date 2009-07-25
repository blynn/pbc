// Type G pairings.

// Requires:
// * mnt.h
// * fops.h
// * param.h
#ifndef __PBC_G_PARAM_H__
#define __PBC_G_PARAM_H__

void pbc_param_init_g(pbc_param_ptr p, fetch_ops_t fops, void *ctx);

/*@manual gparam
Type G curves are generated using the complex multiplication (CM) method.  This
function sets 'p' to a type G pairing parameters from CM parameters 'cm'.
Another part of the library searches for appropriate CM parameters (see below)
and the results can be passed to this function.

To be secure, generic discrete log algorithms must be infeasible in groups of
order r, and finite field discrete log algorithms must be infeasible in finite
fields of order q^6^.  For usual CM parameters, r is a few bits smaller than q.
*/
void pbc_param_init_g_gen(pbc_param_t p, cm_info_ptr cm);

#endif //__PBC_G_PARAM_H__
