// requires
// * stdio.h
// * gmp.h
// * mnt.h
// * fops.h
// * pairing.h
#ifndef __PBC_G_PARAM_H__
#define __PBC_G_PARAM_H__

struct g_param_s {
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

typedef struct g_param_s g_param_t[1];
typedef struct g_param_s *g_param_ptr;

/*@manual gparam
Initialize ''p''. This must be called before ''p'' can be used.
*/
void g_param_init(g_param_ptr p);

/*@manual gparam
Clear ''p''. This should be called after ''p'' is no longer needed.
*/
void g_param_clear(g_param_ptr p);

/*@manual gparam
Write the parameters in ''p'' in a text format onto ''stream''.
*/
void g_param_out_str(FILE *stream, g_param_ptr p);

void g_param_inp_generic (g_param_ptr p, fetch_ops_t fops, void *ctx);
void pairing_init_g_param(pairing_t pairing, g_param_t param);

static inline void g_param_init_inp_str(g_param_ptr p, FILE *stream)
{
    g_param_init(p);
    g_param_init_inp_str(p, stream);
}

/*@manual gparam
Type G curves are generated using the
complex multiplication (CM) method.
This function sets ''p'' to
a type G pairing parameters from CM parameters ''cm''.
Another part of the library searches for
appropriate CM parameters (see below)
and the results can be passed to this function.
</para>
<para>
To be secure, generic discrete log algorithms must
be infeasible in groups of order r, and finite field discrete log algorithms
must be infeasible in finite fields of order q^6.
For usual CM parameters, r is a few bits smaller than q.
*/
void g_param_from_cm(g_param_t p, cm_info_ptr cm);
#endif //__PBC_G_PARAM_H__
