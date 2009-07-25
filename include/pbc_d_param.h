// Type D pairings, aka MNT curves.

// Requires:
// * mnt.h
// * param.h
#ifndef __PBC_D_PARAM_H__
#define __PBC_D_PARAM_H__

void pbc_param_init_d(pbc_param_ptr par, const char *(*tab)(const char *));

/*@manual dparam
Type D curves are generated using the complex multiplication (CM) method.  This
function sets 'p' to a type D pairing parameters from CM parameters 'cm'.
Other library calls search for appropriate CM parameters and the results
can be passed to this function.

To be secure, generic discrete log algorithms must be infeasible in groups of
order r, and finite field discrete log algorithms must be infeasible in finite
fields of order q^6^.  For usual CM parameters, r is a few bits smaller than q.
*/
void pbc_param_init_d_gen(pbc_param_ptr p, cm_info_ptr cm);

#endif //__PBC_D_PARAM_H__
