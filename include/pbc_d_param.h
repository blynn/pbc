// requires
// * stdio.h
// * gmp.h
// * mnt.h
// * fops.h
// * pairing.h
#ifndef __PBC_D_PARAM_H__
#define __PBC_D_PARAM_H__

struct d_param_s {
    mpz_t q; //curve defined over F_q
    mpz_t n; //has order n (= q - t + 1) in F_q
    mpz_t h; //h * r = n, r is prime
    mpz_t r;
    mpz_t a, b; //curve equation is y^2 = x^3 + ax + b
    int k; //embedding degree
    mpz_t nk; //order of curve over F_q^k
    mpz_t hk; //hk * r^2 = nk
    mpz_t *coeff; //coefficients of polynomial used to extend F_q by k/2
    mpz_t nqr; //a quadratic nonresidue in F_q^d that lies in F_q
};

typedef struct d_param_s d_param_t[1];
typedef struct d_param_s *d_param_ptr;

/*@manual dparam
Initialize 'p'. This must be called before 'p' can be used.
*/
void d_param_init(d_param_ptr p);

/*@manual dparam
Clear 'p'. This should be called after 'p' is no longer needed.
*/
void d_param_clear(d_param_ptr p);

/*@manual dparam
Write the parameters in 'p' in a text format onto 'stream'.
*/
void d_param_out_str(FILE *stream, d_param_ptr p);

void d_param_inp_generic (d_param_ptr p, fetch_ops_t fops, void *ctx);
void pairing_init_d_param(pairing_t pairing, d_param_t param);

static inline void d_param_init_inp_str(d_param_ptr p, FILE *stream)
{
    d_param_init(p);
    d_param_init_inp_str(p, stream);
}

/*@manual dparam
Type D curves are generated using the complex multiplication (CM) method.  This
function sets 'p' to a type D pairing parameters from CM parameters 'cm'.
Other library calls search for appropriate CM parameters and the results
can be passed to this function.

To be secure, generic discrete log algorithms must be infeasible in groups of
order r, and finite field discrete log algorithms must be infeasible in finite
fields of order q^6^.  For usual CM parameters, r is a few bits smaller than q.
*/
void d_param_from_cm(d_param_t p, cm_info_ptr cm); #endif //__PBC_D_PARAM_H__
