// Complex floats.
// Called mpc_t, these complex numbers are built on GMP's mpf_t type.

// Requires:
// * stdio.h
// * gmp.h

#ifndef __PBC_MPC_H__
#define __PBC_MPC_H__

#pragma GCC visibility push(hidden)

struct mpc_s {
  mpf_t a;
  mpf_t b;
};
typedef struct mpc_s mpc_t[1];
typedef struct mpc_s *mpc_ptr;

static inline void mpc_init(mpc_ptr c) {
  mpf_init(c->a);
  mpf_init(c->b);
}

static inline void mpc_clear(mpc_ptr c) {
  mpf_clear(c->a);
  mpf_clear(c->b);
}

static inline mpf_ptr mpc_re(mpc_ptr c) {
  return c->a;
}

static inline mpf_ptr mpc_im(mpc_ptr c) {
  return c->b;
}

static inline void mpc_add(mpc_ptr res, mpc_ptr z0, mpc_ptr z1) {
  mpf_add(res->a, z0->a, z1->a);
  mpf_add(res->b, z0->b, z1->b);
}

static inline void mpc_sub(mpc_ptr res, mpc_ptr z0, mpc_ptr z1) {
  mpf_sub(res->a, z0->a, z1->a);
  mpf_sub(res->b, z0->b, z1->b);
}

static inline void mpc_neg(mpc_ptr res, mpc_ptr z) {
  mpf_neg(res->a, z->a);
  mpf_neg(res->b, z->b);
}

static inline void mpc_conj(mpc_ptr res, mpc_ptr z) {
  mpf_set(res->a, z->a);
  mpf_neg(res->b, z->b);
}

static inline void mpc_set(mpc_t res, mpc_t z) {
  mpf_set(res->a, z->a);
  mpf_set(res->b, z->b);
}

static inline void mpc_set_ui(mpc_t res, unsigned long int n) {
  mpf_set_ui(res->a, n);
  mpf_set_ui(res->b, 0);
}

static inline void mpc_add_ui(mpc_t res, mpc_t z, unsigned long int n) {
  mpf_add_ui(res->a, z->a, n);
}

static inline void mpc_mul_ui(mpc_t res, mpc_t z, unsigned long int n) {
  mpf_mul_ui(res->a, z->a, n);
  mpf_mul_ui(res->b, z->b, n);
}

static inline void mpc_mul_mpf(mpc_t res, mpc_t z, mpf_t f) {
  mpf_mul(res->a, z->a, f);
  mpf_mul(res->b, z->b, f);
}

void mpc_mul(mpc_t res, mpc_t z0, mpc_t z1);
void mpc_mul_2exp(mpc_t res, mpc_t z, unsigned long int);
void mpc_div(mpc_t res, mpc_t z0, mpc_t z1);
void mpc_muli(mpc_t res, mpc_t z);
void mpc_sqr(mpc_t res, mpc_t z);
void mpc_inv(mpc_t res, mpc_t z);
size_t mpc_out_str(FILE *stream, int base, size_t n_digits, mpc_t op);
void mpc_pow_ui(mpc_t res, mpc_t z, unsigned int n);

#pragma GCC visibility pop

#endif //__PBC_MPC_H__
