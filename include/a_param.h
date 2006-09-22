#ifndef A_PARAM_H
#define A_PARAM_H

#include "fieldquadratic.h"
#include "pairing.h"
#include "fops.h"

struct a_param_s {
    int exp2;
    int exp1;
    int sign1;
    int sign0;
    mpz_t r; // r = 2^exp2 + sign1 * 2^exp1 + sign0 * 1
    mpz_t q; // we work in E(F_q) (and E(F_q^2))
    mpz_t h; // r * h = q + 1
};
typedef struct a_param_s a_param_t[1];
typedef struct a_param_s *a_param_ptr;

void a_param_init(a_param_t p);
void a_param_clear(a_param_t p);
void a_param_gen(a_param_t p, int rbits, int qbits);
void a_param_out_str(FILE *stream, a_param_ptr p);
void a_param_inp_generic (a_param_ptr p, fetch_ops_t *fops, void *ctx);
void a_param_inp_buf(a_param_ptr p, const char *buf, size_t len);
void a_param_inp_str(a_param_ptr p, FILE *stream);
void pairing_init_a_param(pairing_t pairing, a_param_t p);

#endif //A_PARAM_H
