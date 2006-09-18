#ifndef BGN_PARAM_H
#define BGN_PARAM_H

#include "fieldquadratic.h"
#include "pairing.h"
#include "fops.h"

struct bgn_param_s {
    mpz_t p;
    mpz_t n;
    unsigned int l;
};
typedef struct bgn_param_s bgn_param_t[1];
typedef struct bgn_param_s *bgn_param_ptr;

void bgn_param_init(bgn_param_t param);
void bgn_param_clear(bgn_param_t param);
void bgn_param_out_str(FILE *stream, bgn_param_ptr p);
void bgn_param_inp_generic(bgn_param_ptr p, fetch_ops_t *fops, void *ctx);
void bgn_param_inp_buf(bgn_param_ptr p, const char *buf, size_t len);
void bgn_param_inp_str(bgn_param_ptr p, FILE *stream);
void bgn_param_gen(bgn_param_t param, mpz_t order);
void pairing_init_bgn_param(pairing_t pairing, bgn_param_t param);

#endif //BGN_PARAM_H
