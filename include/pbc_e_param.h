// Type E pairings.

// Requires:
// * param.h
#ifndef __PBC_E_PARAM_H__
#define __PBC_E_PARAM_H__

struct symtab_s;
int pbc_param_init_e(pbc_param_ptr par, struct symtab_s *tab);

/*@manual eparam
Generate type E pairing parameters and store them in 'p',
where the group order r is 'rbits' long, and the order of the base field q
is 'qbits' long. To be secure, generic discrete log algorithms must
be infeasible in groups of order r, and finite field discrete log algorithms
must be infeasible in finite fields of order q,
e.g. 'rbits' = 160, 'qbits' = 1024.

This pairing is just a curiosity: it can be implemented entirely in a field of
prime order, that is, only arithmetic modulo a prime is needed and there is
never a need to extend a field.

If discrete log in field extensions are found to be substantially easier to
solve than previously thought, or discrete log can be solved in elliptic curves
as easily as they can be in finite fields, this pairing type may become useful.
*/
void pbc_param_init_e_gen(pbc_param_t p, int rbits, int qbits);

#endif //__PBC_E_PARAM_H__
