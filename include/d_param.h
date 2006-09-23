#ifndef C_PARAM_H
#define C_PARAM_H

#include "pairing.h"
#include "mnt.h"
#include "fops.h"

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

void d_param_init(d_param_ptr cc);
void d_param_clear(d_param_ptr cc);
void d_param_out_str(FILE *stream, d_param_ptr p);
void d_param_inp_generic (d_param_ptr p, fetch_ops_t fops, void *ctx);
void d_param_inp_buf (d_param_ptr p, const char *buf, size_t len);
void d_param_inp_str(d_param_ptr p, FILE *stream);
void pairing_init_c_param(pairing_t pairing, d_param_t param);

static inline void d_param_init_inp_str(d_param_ptr p, FILE *stream)
{
    d_param_init(p);
    d_param_init_inp_str(p, stream);
}

void d_param_from_cm(d_param_t param, cm_info_ptr cm);
#endif //C_PARAM_H
