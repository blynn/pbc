//pairing-based signatures library
#include <string.h> //for memcmp()

#include "pbc_sig.h"
#include "hash.h"

void bls_sign(unsigned char *sig, unsigned int hashlen, unsigned char *hash,
	bls_private_key_t sk)
{
    element_t h;

    element_init_G1(h, sk->param->pairing);
    element_from_hash(h, hashlen, hash);
    element_pow_zn(h, h, sk->x);
    element_to_bytes_x_only(sig, h);

    element_clear(h);
}

int bls_verify(unsigned char *sig, unsigned int hashlen, unsigned char *hash,
	bls_public_key_t pk)
{
    //have to mess with internals since we are only given the x-coord
    element_t hx;
    element_t h;
    int res;

    pairing_ptr pairing = pk->param->pairing;

    element_init_G1(h, pairing);
    element_from_hash(h, hashlen, hash);
    element_init_G1(hx, pairing);
    element_from_bytes_x_only(hx, sig);

    res = is_almost_coddh(h, hx, pk->param->g, pk->gx, pk->param->pairing);
    element_clear(hx);
    return res;
}

void bls_gen_sys_param(bls_sys_param_t param, pairing_t pairing)
{
    param->pairing = pairing;
    element_init_G2(param->g, pairing);
    element_random(param->g);
    param->signature_length = pairing_length_in_bytes_x_only_G1(pairing);
}

void bls_gen(bls_public_key_t pk, bls_private_key_t sk, bls_sys_param_t param)
{
    pk->param = sk->param = param;
    element_init_G2(pk->gx, param->pairing);
    element_init_Zr(sk->x, param->pairing);
    element_random(sk->x);
    element_pow_zn(pk->gx, param->g, sk->x);
}

void bb_gen_sys_param(bb_sys_param_t param, pairing_t pairing)
{
    param->pairing = pairing;
    //signature is a point (only need one coordinate) and an element of Z_r
    param->signature_length = pairing_length_in_bytes_x_only_G1(pairing)
	    + pairing_length_in_bytes_Zr(pairing);
}

void bb_gen(bb_public_key_t pk, bb_private_key_t sk, bb_sys_param_t param)
{
    pairing_ptr pairing = param->pairing;
    pk->param = sk->param = param;
    
    element_init(sk->x, pairing->Zr);
    element_init(sk->y, pairing->Zr);
    element_random(sk->x);
    element_random(sk->y);
    element_init(pk->g1, param->pairing->G1);
    element_init(pk->g2, param->pairing->G2);
    element_init(pk->z, param->pairing->GT);
    element_random(pk->g2);
    element_random(pk->g1);
    element_init(pk->u, param->pairing->G2);
    element_init(pk->v, param->pairing->G2);
    element_pow_zn(pk->u, pk->g2, sk->x);
    element_pow_zn(pk->v, pk->g2, sk->y);
    bilinear_map(pk->z, pk->g1, pk->g2, param->pairing);
}

void bb_sign(unsigned char *sig, unsigned int hashlen, unsigned char *hash,
	bb_public_key_t pk, bb_private_key_t sk)
{
    int len;
    element_t sigma;
    element_t r, z, m;
    bb_sys_param_ptr param = pk->param;
    pairing_ptr pairing = param->pairing;

    element_init(r, pairing->Zr);
    element_init(z, pairing->Zr);
    element_init(m, pairing->Zr);

    element_random(r);
    mpz_import(m->data, hashlen, 1, 1, 1, 0, hash);
    mpz_mod(m->data, m->data, pairing->r);
    element_mul(z, sk->y, r);
    element_add(z, z, sk->x);
    element_add(z, z, m);
    element_invert(z, z);
    element_init(sigma, pairing->G1);
    element_pow_zn(sigma, pk->g1, z);

    len = element_to_bytes_x_only(sig, sigma);
    element_to_bytes(&sig[len], r);

    element_clear(sigma);
    element_clear(r);
    element_clear(z);
    element_clear(m);
}

int bb_verify(unsigned char *sig, unsigned int hashlen, unsigned char *hash,
	bb_public_key_t pk)
{
    element_t sigma, r;
    element_t m;
    element_t t0, t1, t2;
    int res;
    int len;
    pairing_ptr pairing = pk->param->pairing;

    element_init(m, pairing->Zr);
    mpz_import(m->data, hashlen, 1, 1, 1, 0, hash);
    mpz_mod(m->data, m->data, pairing->r);

    element_init(sigma, pairing->G1);
    len = element_from_bytes_x_only(sigma, sig);

    element_init(r, pairing->Zr);
    element_from_bytes(r, sig + len);

    element_init(t0, pairing->G2);
    element_init(t1, pairing->G2);
    element_init(t2, pairing->GT);

    element_pow_zn(t0, pk->g2, m);
    element_pow_zn(t1, pk->v, r);
    element_mul(t0, t0, t1);
    element_mul(t0, t0, pk->u);
    bilinear_map(t2, sigma, t0, pairing);
    if (!element_cmp(t2, pk->z)) {
	res = 1;
    } else {
	element_mul(t2, t2, pk->z);
	res = element_is1(t2);
    }

    element_clear(t0);
    element_clear(t1);
    element_clear(t2);
    element_clear(m);
    element_clear(sigma);
    element_clear(r);
    return res;
}

void ib_setup(ib_sys_param_t param, ib_master_key_t master, pairing_t pairing)
{
    param->pairing = pairing;
    element_init(param->g, pairing->G2);
    element_init(param->gx, pairing->G2);
    element_random(param->g);
    element_init(master->x, pairing->Zr);
    element_random(master->x);
    element_pow_zn(param->gx, param->g, master->x);
    master->param = param;
}

void ib_extract(ib_private_key_t priv, unsigned int idlen, unsigned char *id,
	ib_master_key_t sk)
{
    priv->param = sk->param;

    element_init(priv->d, sk->param->pairing->G1);
    element_init(priv->q, sk->param->pairing->G1);
    element_from_hash(priv->q, idlen, id);
    element_pow_zn(priv->d, priv->q, sk->x);
}

int ib_keyagree_length(ib_sys_param_t param)
{
    return param->pairing->GT->fixed_length_in_bytes;
}

void ib_keyagree(unsigned char *key, ib_private_key_t priv, unsigned int idlen, unsigned char *id)
    //requires G1 == G2
{
    element_t u;
    element_t e;

    element_init(u, priv->param->pairing->G1);
    element_init(e, priv->param->pairing->GT);
    element_from_hash(u, idlen, id);
    bilinear_map(e, u, priv->d, priv->param->pairing);

    element_to_bytes(key, e);
    element_clear(e);
    element_clear(u);
}

//TODO: only use x-coord
void cc_sign(unsigned char *sig, unsigned int hashlen, unsigned char *hash,
	ib_private_key_t sk)
{
    unsigned char *buf, digest[hash_length];
    element_t r, h;
    element_t u, v;
    hash_ctx_t ctx;
    unsigned int len;

    element_init(h, sk->param->pairing->Zr);

    element_init(r, sk->param->pairing->Zr);
    element_random(r);

    element_init(v, sk->param->pairing->G1);

    element_init(u, sk->param->pairing->G1);
    element_pow_zn(u, sk->q, r);

    buf = malloc(len = element_length_in_bytes(u));
    element_to_bytes(buf, u);
    hash_init(ctx);
    hash_update(ctx, buf, len);
    hash_update(ctx, hash, hashlen);
    hash_final(digest, ctx);
    free(buf);

    element_from_hash(h, hash_length, digest);
    element_add(r, r, h);
    element_pow_zn(v, sk->d, r);
    len = element_to_bytes(sig, u);
    element_to_bytes(&sig[len], v);

    element_clear(r);
    element_clear(h);
    element_clear(u);
    element_clear(v);
}

int cc_signature_length(ib_sys_param_t param)
{
    return 2 * param->pairing->G1->fixed_length_in_bytes;
}

int cc_verify(unsigned char *sig, unsigned int hashlen, unsigned char *hash,
	unsigned int idlen, unsigned char *id, ib_sys_param_t param)
{
    int len;
    int res;
    element_t u, v, temp;
    element_t h;
    unsigned char digest[hash_length];
    hash_ctx_t ctx;

    element_init(u, param->pairing->G1);
    element_init(v, param->pairing->G1);
    element_init(temp, param->pairing->G1);

    len = element_from_bytes(u, sig);
    element_from_bytes(v, &sig[len]);

    hash_init(ctx);
    hash_update(ctx, sig, len);
    hash_update(ctx, hash, hashlen);
    hash_final(digest, ctx);

    element_init(h, param->pairing->Zr);

    element_from_hash(h, hash_length, digest);
    element_from_hash(temp, idlen, id);
    element_pow_zn(temp, temp, h);
    element_mul(temp, temp, u);

    res = is_almost_coddh(temp, v, param->g, param->gx, param->pairing);

    element_clear(u);
    element_clear(v);
    element_clear(temp);
    element_clear(h);
    return res;
}

// This scheme is patented.
/*
//TODO: only use x-coord
void skschnorr_sign(unsigned char *sig, unsigned int hashlen, unsigned char *hash,
	ib_private_key_t sk)
{
    unsigned char *buf;
    element_t h;
    element_t r;
    element_t u, temp, v;
    hash_ctx_t ctx;
    unsigned int len;

    element_t e;
    element_init(u, sk->param->pairing->G1);
    element_init(temp, sk->param->pairing->G1);
    element_init(e, sk->param->pairing->GT);
    element_init(h, sk->param->pairing->Zr);
    element_init(r, sk->param->pairing->Zr);
    element_random(r);

    element_init(v, sk->param->pairing->G2);
    element_pow_zn(v, sk->param->g, r);
    bilinear_map(e, sk->q, v, sk->param->pairing);

    buf = malloc(len = element_length_in_bytes(e));
    element_to_bytes(buf, e);
    hash_init(ctx);
    hash_update(ctx, hash, hashlen);
    hash_update(ctx, buf, len);
    hash_final(sig, ctx);
    free(buf);

    element_from_hash(h, hash_length, sig);
    element_pow_zn(u, sk->d, h);
    element_pow_zn(temp, sk->q, r);
    element_mul(u, u, temp);
    element_to_bytes(&sig[hash_length], u);

    element_clear(h);
    element_clear(r);
    element_clear(u);
    element_clear(temp);
    element_clear(v);
}

int skschnorr_signature_length(ib_sys_param_t param)
{
    return param->pairing->G1->fixed_length_in_bytes + hash_length;
}

int skschnorr_verify(unsigned char *sig, unsigned int hashlen, unsigned char *hash,
	unsigned int idlen, unsigned char *id, ib_sys_param_t param)
{
    int res;
    element_t e, e2;
    element_t u, v;
    element_t h;
    unsigned char digest[hash_length];
    hash_ctx_t ctx;
    unsigned int len;
    unsigned char *buf;

    element_init(e, param->pairing->GT);
    element_init(e2, param->pairing->GT);
    element_init(h, param->pairing->Zr);
    element_init(u, param->pairing->G1);
    element_init(v, param->pairing->G2);

    element_from_hash(h, hash_length, sig);
    element_neg(h, h);

    element_pow_zn(v, param->gx, h);

    element_from_bytes(u, &sig[hash_length]);
    bilinear_map(e, u, param->g, param->pairing);

    element_from_hash(u, idlen, id);
    bilinear_map(e2, u, v, param->pairing);

    element_mul(e, e, e2);

    buf = malloc(len = element_length_in_bytes(e));
    element_to_bytes(buf, e);
    hash_init(ctx);
    hash_update(ctx, hash, hashlen);
    hash_update(ctx, buf, len);
    hash_final(digest, ctx);

    res = !memcmp(digest, sig, hash_length);

    element_clear(u);
    element_clear(v);
    element_clear(e);
    element_clear(e2);
    element_clear(h);
    return res;
}
*/
