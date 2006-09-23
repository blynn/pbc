#ifndef F_PARAM_H
#define F_PARAM_H

#include "fops.h"

struct f_param_s {
    mpz_t q; //curve defined over F_q
    mpz_t r; //r is the order of the curve
    mpz_t b; //curve equation is y^2 = x^3 + b
    mpz_t beta; //beta is a quadratic nonresidue in Fq
	//we use F_q^2 = F_q[sqrt(beta)]
    mpz_t alpha0, alpha1;
	//the polynomial x^6 + alpha0 + alpha1 sqrt(beta)
	//is irreducible over F_q^2[x], so
	//we can extend F_q^2 to F_q^12 using the
	//sixth root of -(alpha0 + alpha1 sqrt(beta))
};
typedef struct f_param_s f_param_t[1];
typedef struct f_param_s *f_param_ptr;

void f_param_init(f_param_t fp);
void f_param_clear(f_param_t fp);
void f_param_gen(f_param_t fp, int bits);
void f_param_out_str(FILE *stream, f_param_ptr p);
void f_param_inp_generic (f_param_ptr p, fetch_ops_t fops, void *ctx);
void f_param_inp_buf (f_param_ptr p, const char *buf, size_t len);
void f_param_inp_str (f_param_ptr p, FILE *stream);
void pairing_init_f_param(pairing_t pairing, f_param_t param);

#endif //F_PARAM_H
