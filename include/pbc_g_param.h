// Type G pairings.

// Requires:
// * mnt.h
// * param.h
#ifndef __PBC_G_PARAM_H__
#define __PBC_G_PARAM_H__

struct symtab_s;
int pbc_param_init_g(pbc_param_ptr par, struct symtab_s *tab);

/*@manual gparam
Type G curves are generated using the complex multiplication (CM) method.  This
function sets 'p' to a type G pairing parameters from CM parameters 'cm'.
They have embedding degree 10.

To be secure, generic discrete log algorithms must be infeasible in groups of
order r, and finite field discrete log algorithms must be infeasible in finite
fields of order q^6^.  For usual CM parameters, r is a few bits smaller than q.

They are quite slow at the moment so for now type F is a better choice.

The file `param/g149.param` contains parameters for a
type G pairing with 149-bit group and field sizes.
*/
void pbc_param_init_g_gen(pbc_param_t p, pbc_cm_ptr cm);

#endif //__PBC_G_PARAM_H__
