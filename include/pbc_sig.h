#ifndef SIG_H
#define SIG_H

//pairing-based signatures library

#include "pbc.h"

struct bls_sys_param_s {
    pairing_ptr pairing;
    element_t g;
    int signature_length;
};
typedef struct bls_sys_param_s bls_sys_param_t[1];
typedef struct bls_sys_param_s *bls_sys_param_ptr;

struct bls_private_key_s {
    bls_sys_param_ptr param;
    element_t x;
};
typedef struct bls_private_key_s bls_private_key_t[1];
typedef struct bls_private_key_s *bls_private_key_ptr;

struct bls_public_key_s {
    bls_sys_param_ptr param;
    element_t gx;
};
typedef struct bls_public_key_s bls_public_key_t[1];
typedef struct bls_public_key_s *bls_public_key_ptr;

void bls_gen_sys_param(bls_sys_param_t param, pairing_t pairing);
void bls_gen(bls_public_key_t pk, bls_private_key_t sk, bls_sys_param_t param);
void bls_sign(unsigned char *sig, unsigned int hashlen, unsigned char *hash,
	bls_private_key_t sk);
int bls_verify(unsigned char *sig, unsigned int hashlen, unsigned char *hash,
	bls_public_key_t pk);

struct bb_sys_param_s {
    pairing_ptr pairing;
    int signature_length;
};
typedef struct bb_sys_param_s bb_sys_param_t[1];
typedef struct bb_sys_param_s *bb_sys_param_ptr;

struct bb_private_key_s {
    bb_sys_param_ptr param;
    element_t x, y;
};
typedef struct bb_private_key_s bb_private_key_t[1];
typedef struct bb_private_key_s *bb_private_key_ptr;

struct bb_public_key_s {
    bb_sys_param_ptr param;
    element_t g1, g2, u, v, z;
};
typedef struct bb_public_key_s bb_public_key_t[1];
typedef struct bb_public_key_s *bb_public_key_ptr;

void bb_gen_sys_param(bb_sys_param_t param, pairing_t pairing);
void bb_gen(bb_public_key_t pk, bb_private_key_t sk, bb_sys_param_t param);
void bb_sign(unsigned char *sig, unsigned int hashlen, unsigned char *hash,
	bb_public_key_t pk, bb_private_key_t sk);
int bb_verify(unsigned char *sig, unsigned int hashlen, unsigned char *hash,
	bb_public_key_t pk);

struct ib_sys_param_s {
    pairing_ptr pairing;
    element_t g, gx;
};
typedef struct ib_sys_param_s ib_sys_param_t[1];
typedef struct ib_sys_param_s *ib_sys_param_ptr;

struct ib_master_key_s {
    ib_sys_param_ptr param;
    element_t x;
};
typedef struct ib_master_key_s ib_master_key_t[1];
typedef struct ib_master_key_s *ib_master_key_ptr;

struct ib_private_key_s {
    ib_sys_param_ptr param;
    element_t q;
    element_t d;
};
typedef struct ib_private_key_s ib_private_key_t[1];
typedef struct ib_private_key_s *ib_private_key_ptr;

void ib_setup(ib_sys_param_t param, ib_master_key_t pkgpriv, pairing_t pairing);
void ib_extract(ib_private_key_t priv, unsigned int idlen, unsigned char *id,
	ib_master_key_t sk);
int cc_signature_length(ib_sys_param_t param);
void cc_sign(unsigned char *sig, unsigned int hashlen, unsigned char *hash,
	ib_private_key_t sk);
int cc_verify(unsigned char *sig, unsigned int hashlen, unsigned char *hash,
	unsigned int idlen, unsigned char *id, ib_sys_param_t param);

/* This scheme is patented
void skschnorr_sign(unsigned char *sig, unsigned int hashlen, unsigned char *hash,
	ib_private_key_t sk);
int skschnorr_signature_length(ib_sys_param_t param);
int skschnorr_verify(unsigned char *sig, unsigned int hashlen, unsigned char *hash,
	unsigned int idlen, unsigned char *id, ib_sys_param_t param);
*/

#endif //SIG_H
