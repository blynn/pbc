// requires
// * stdio.h
// * gmp.h
// * pairing.h
// * fops.h
#ifndef __PBC_E_PARAM_H__
#define __PBC_E_PARAM_H__

struct e_param_s {
    mpz_t q; //curve defined over F_q
    mpz_t r; //q = h r^2 + 1, r is prime
    mpz_t h; //h is 28 h'^2 for some h'
    mpz_t a, b; //curve equation is y^2 = x^3 + ax + b
    int exp2;
    int exp1;
    int sign1;
    int sign0;
};
typedef struct e_param_s e_param_t[1];
typedef struct e_param_s *e_param_ptr;

/*@manual eparam
Initialize ''p''. This must be called before ''p'' can be used.
*/
void e_param_init(e_param_t ep);

/*@manual eparam
Clear ''p''. This should be called after ''p'' is no longer needed.
*/
void e_param_clear(e_param_t ep);

/*@manual eparam
Generate type E pairing parameters and store them in ''p'',
where the group order r is ''rbits'' long, and the order of the base field q
is ''qbits'' long. To be secure, generic discrete log algorithms must
be infeasible in groups of order r, and finite field discrete log algorithms
must be infeasible in finite fields of order q.
Typical values: ''rbits'' = 160, ''qbits'' = 1024.
*/
void e_param_gen(e_param_t p, int rbits, int qbits);

/*@manual eparam
Write the parameters in ''p'' in a text format onto ''stream''.
*/
void e_param_out_str(FILE *stream, e_param_ptr p);

void e_param_inp_generic (e_param_ptr p, fetch_ops_t fops, void *ctx);
void pairing_init_e_param(pairing_t pairing, e_param_t param);

#endif //__PBC_E_PARAM_H__
