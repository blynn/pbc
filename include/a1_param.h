#ifndef A1_PARAM_H
#define A1_PARAM_H

#include "fieldquadratic.h"
#include "pairing.h"
#include "fops.h"

struct a1_param_s {
    mpz_t p;
    mpz_t n;
    unsigned int l;
};
typedef struct a1_param_s a1_param_t[1];
typedef struct a1_param_s *a1_param_ptr;

void a1_param_init(a1_param_t param);
void a1_param_clear(a1_param_t param);
void a1_param_out_str(FILE *stream, a1_param_ptr p);
void a1_param_inp_generic(a1_param_ptr p, fetch_ops_t fops, void *ctx);
void a1_param_inp_buf(a1_param_ptr p, const char *buf, size_t len);
void a1_param_inp_str(a1_param_ptr p, FILE *stream);
void a1_param_gen(a1_param_t param, mpz_t order);
void pairing_init_a1_param(pairing_t pairing, a1_param_t param);

#endif //A1_PARAM_H
