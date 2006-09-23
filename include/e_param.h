#ifndef E_PARAM_H
#define E_PARAM_H

#include "pairing.h"
#include "fops.h"

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

void e_param_init(e_param_t ep);
void e_param_clear(e_param_t ep);
void e_param_gen(e_param_t p, int rbits, int qbits);
void e_param_out_str(FILE *stream, e_param_ptr p);
void e_param_inp_generic (e_param_ptr p, fetch_ops_t fops, void *ctx);
void e_param_inp_buf (e_param_ptr p, const char *buf, size_t len);
void e_param_inp_str(e_param_ptr p, FILE *stream);
void pairing_init_e_param(pairing_t pairing, e_param_t param);

#endif //E_PARAM_H
