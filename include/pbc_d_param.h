// Type D pairings, aka MNT curves.

// Requires:
// * mnt.h
// * param.h
#ifndef __PBC_D_PARAM_H__
#define __PBC_D_PARAM_H__

struct symtab_s;
int pbc_param_init_d(pbc_param_ptr par, struct symtab_s *tab);

/*@manual dparam
Type D curves are generated using the complex multiplication (CM) method.  This
function sets 'p' to a type D pairing parameters from CM parameters 'cm'.
Other library calls search for appropriate CM parameters and the results
can be passed to this function.

To be secure, generic discrete log algorithms must be infeasible in groups of
order r, and finite field discrete log algorithms must be infeasible in finite
fields of order q^6^.  For usual CM parameters, r is a few bits smaller than q.

Using type D pairings allows elements of group G1 to be quite short, typically
170-bits. Because of a certain trick, elements of group G2 need only be 3 times
longer, that is, about 510 bits rather than 6 times long. They are not quite
as short as type F pairings, but much faster.

I sometimes refer to a type D curve as a triplet of numbers: the discriminant,
the number of bits in the prime q, and the number of bits in the prime r. The
`gen/listmnt` program prints these numbers.

Among the bundled type D curve parameters are the curves 9563-201-181,
62003-159-158 and 496659-224-224 which have shortened names `param/d201.param`,
`param/d159.param` and `param/d225.param` respectively.

See `gen/listmnt.c` and `gen/gendparam.c` for how to generate type D pairing
parameters.
*/
void pbc_param_init_d_gen(pbc_param_ptr p, pbc_cm_ptr cm);

#endif //__PBC_D_PARAM_H__
