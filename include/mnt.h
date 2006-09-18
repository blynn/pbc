#ifndef MNT_H
#define MNT_H

#include "darray.h"
#include "curve.h"

struct cm_info_s {
    mpz_t q; //curve defined over F_q
    mpz_t n; //has order n (= q - t + 1) in F_q (and r^2 in F_q^k)
    mpz_t h; //h * r = n, r is prime
    mpz_t r;
    int D; //discrminant needed to find j-invariant
    int k; //embedding degree
};

typedef struct cm_info_s *cm_info_ptr;
typedef struct cm_info_s cm_info_t[1];

void cm_info_init(cm_info_t cm);
void cm_info_clear(cm_info_t cm);

int find_mnt6_curve(darray_t L, unsigned int D, unsigned int bitlimit);

#endif //MNT_H
