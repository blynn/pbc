//To understand this program,
//see Boneh, Boyen and Shacham, "Short Group Signatures"
#include "pbc.h"
#include "get_time.h"
#include "utils.h"

struct bbs_sys_param_s {
    pairing_ptr pairing;
    int signature_length;
};
typedef struct bbs_sys_param_s bbs_sys_param_t[1];
typedef struct bbs_sys_param_s *bbs_sys_param_ptr;

struct bbs_group_public_key_s {
    bbs_sys_param_ptr param;
    element_t g1, g2;
    element_t h, u, v, w;
    /* and precomputed values */
    element_t pr_g1_g2;
    element_t pr_h_g2;
    element_t pr_h_w;
    element_t pr_g1_g2_inv;
};
typedef struct bbs_group_public_key_s bbs_group_public_key_t[1];
typedef struct bbs_group_public_key_s *bbs_group_public_key_ptr;

struct bbs_group_private_key_s {
    bbs_sys_param_ptr param;
    element_t A;
    element_t x;
    /* and precomputed values */
    element_t pr_A_g2;
};
typedef struct bbs_group_private_key_s bbs_group_private_key_t[1];
typedef struct bbs_group_private_key_s *bbs_group_private_key_ptr;

struct bbs_manager_private_key_s {
    bbs_sys_param_ptr param;
    mpz_t xi1, xi2;
};
typedef struct bbs_manager_private_key_s bbs_manager_private_key_t[1];
typedef struct bbs_manager_private_key_s *bbs_manager_private_key_ptr;

void bbs_gen_sys_param(bbs_sys_param_t param, pairing_t pairing)
{
    param->pairing = pairing;
    param->signature_length = 3 * pairing->G1->fixed_length_in_bytes
	+ 6 * pairing->Zr->fixed_length_in_bytes;
}

void bbs_gen(bbs_group_public_key_t gpk, bbs_manager_private_key_t gmsk,
	int n, bbs_group_private_key_t *gsk, bbs_sys_param_t param)
{
    pairing_ptr pairing = param->pairing;
    mpz_t z0;
    mpz_t gamma;
    mpz_ptr r = pairing->r;
    int i;

    gpk->param = param;
    gmsk->param = param;
    element_init(gpk->g1, pairing->G1);
    element_init(gpk->g2, pairing->G2);
    element_init(gpk->h, pairing->G1);
    element_init(gpk->u, pairing->G1);
    element_init(gpk->v, pairing->G1);
    element_init(gpk->w, pairing->G2);
    mpz_init(gmsk->xi1);
    mpz_init(gmsk->xi2);
    mpz_init(gamma);
    mpz_init(z0);

    element_random(gpk->g2);
    element_random(gpk->g1);
    element_random(gpk->h);
    pbc_mpz_random(gmsk->xi1, r);
    pbc_mpz_random(gmsk->xi2, r);
    mpz_invert(z0, gmsk->xi1, r);
    element_pow(gpk->u, gpk->h, z0);
    mpz_invert(z0, gmsk->xi2, r);
    element_pow(gpk->v, gpk->h, z0);
    pbc_mpz_random(gamma, r);
    element_pow(gpk->w, gpk->g2, gamma);

    for (i=0; i<n; i++) {
	gsk[i]->param = param;
	element_init(gsk[i]->A, pairing->G1);
	element_init(gsk[i]->x, pairing->Zr);

	element_random(gsk[i]->x);
	//TODO: "->data" is bad
	mpz_add(z0, gamma, gsk[i]->x->data);
	mpz_invert(z0, z0, r);
	element_pow(gsk[i]->A, gpk->g1, z0);

        /* do some precomputation */
        /* TODO: could instead compute from e(g1,g2) ... */
        element_init(gsk[i]->pr_A_g2, pairing->GT);
        bilinear_map(gsk[i]->pr_A_g2, gsk[i]->A, gpk->g2, pairing);
    }


    /* do some precomputation */
    element_init(gpk->pr_g1_g2, pairing->GT);
    element_init(gpk->pr_g1_g2_inv, pairing->GT);
    element_init(gpk->pr_h_g2, pairing->GT);
    element_init(gpk->pr_h_w, pairing->GT);
    bilinear_map(gpk->pr_g1_g2, gpk->g1, gpk->g2, pairing);
    element_invert(gpk->pr_g1_g2_inv, gpk->pr_g1_g2);
    bilinear_map(gpk->pr_h_g2, gpk->h, gpk->g2, pairing);
    bilinear_map(gpk->pr_h_w, gpk->h, gpk->w, pairing);


    mpz_clear(z0);
    mpz_clear(gamma);
}

static mpz_ptr to_mpz(element_t e) {
    return e->data;
}

void bbs_sign(unsigned char *sig,
	int hashlen, void *hash,
	bbs_group_public_key_t gpk, bbs_group_private_key_t gsk)
{
    bbs_sys_param_ptr param = gpk->param;
    pairing_ptr pairing = param->pairing;
    field_ptr Fp = pairing->Zr;
    element_t T1, T2, T3;
    element_t R1, R2, R3, R4, R5;
    element_t alpha, beta;
    element_t c;
    element_t ralpha, rbeta, rx, rdelta1, rdelta2;
    element_t z0, z1;
    element_t e10, et0;
    unsigned char *writeptr = sig;
    UNUSED_VAR (hashlen);
    UNUSED_VAR (hash);

    element_init(T1, pairing->G1);
    element_init(T2, pairing->G1);
    element_init(T3, pairing->G1);
    element_init(R1, pairing->G1);
    element_init(R2, pairing->G1);
    element_init(R3, pairing->GT);
    element_init(R4, pairing->G1);
    element_init(R5, pairing->G1);

    element_init(c, Fp);
    element_init(alpha, Fp); element_random(alpha);
    element_init(beta, Fp); element_random(beta);

    //temp variables
    element_init(z0, Fp);
    element_init(z1, Fp);
    element_init(et0, pairing->GT);
    element_init(e10, pairing->G1);

    element_init(ralpha, Fp); element_random(ralpha);
    element_init(rbeta, Fp); element_random(rbeta);
    element_init(rx, Fp); element_random(rx);
    element_init(rdelta1, Fp); element_random(rdelta1);
    element_init(rdelta2, Fp); element_random(rdelta2);

    element_pow(T1, gpk->u, to_mpz(alpha));
    element_pow(T2, gpk->v, to_mpz(beta));
    element_add(z0, alpha, beta);

    element_pow(T3, gpk->h, to_mpz(z0));
    element_mul(T3, T3, gsk->A);

    element_pow(R1, gpk->u, to_mpz(ralpha));

    element_pow(R2, gpk->v, to_mpz(rbeta));

    /*
     * rather than computing e(T3,g2), note that T3 = A h^{alpha+beta},
     * use precomputed e(A,g2) and e(h,g2), and use appropriate
     * exponentiations in GT.
     */

    //bilinear_map(et0, T3, gpk->g2, pairing);  /* precomputed */
    element_pow(et0, gpk->pr_h_g2, to_mpz(z0)); /* NB. here z0 = alpha+beta */
    element_mul(et0, et0, gsk->pr_A_g2);
    //element_pow(R3, et0, to_mpz(rx));

    // bilinear_map(et0, gpk->h, gpk->w, pairing);  /* precomputed */
    element_add(z0, ralpha, rbeta);
    element_neg(z0, z0);
    //element_pow(et0, gpk->pr_h_w, to_mpz(z0));
    //element_mul(R3, R3, et0);
    // bilinear_map(et0, gpk->h, gpk->g2, pairing);  /* precomputed */
    element_add(z1, rdelta1, rdelta2);
    element_neg(z1, z1);
    //element_pow(et0, gpk->pr_h_g2, to_mpz(z1));
    //element_mul(R3, R3, et0);

    element_pow3(R3, et0, to_mpz(rx),
                 gpk->pr_h_w, to_mpz(z0), gpk->pr_h_g2, to_mpz(z1));

    //element_pow(R4, T1, to_mpz(rx));
    element_neg(z0, rdelta1);
    //element_pow(e10, gpk->u, to_mpz(z0));
    //element_mul(R4, R4, e10);
    element_pow2(R4, T1, to_mpz(rx), gpk->u, to_mpz(z0));

    //element_pow(R5, T2, to_mpz(rx));
    element_neg(z0, rdelta2);
    //element_pow(e10, gpk->v, to_mpz(z0));
    //element_mul(R5, R5, e10);
    element_pow2(R5, T2, to_mpz(rx), gpk->v, to_mpz(z0));

    //TODO: c should be the hash of T's and R's
    element_random(c);

    //now the r's represent the values of the s's
    //no need to allocate yet more variables
    element_mul(z0, c, alpha);
    element_add(ralpha, ralpha, z0);

    element_mul(z0, c, beta);
    element_add(rbeta, rbeta, z0);

    element_mul(z1, c, gsk->x);
    element_add(rx, rx, z1);

    element_mul(z0, z1, alpha);
    element_add(rdelta1, rdelta1, z0);

    element_mul(z0, z1, beta);
    element_add(rdelta2, rdelta2, z0);

    writeptr += element_to_bytes(writeptr, T1);
    writeptr += element_to_bytes(writeptr, T2);
    writeptr += element_to_bytes(writeptr, T3);
    writeptr += element_to_bytes(writeptr, c);
    writeptr += element_to_bytes(writeptr, ralpha);
    writeptr += element_to_bytes(writeptr, rbeta);
    writeptr += element_to_bytes(writeptr, rx);
    writeptr += element_to_bytes(writeptr, rdelta1);
    writeptr += element_to_bytes(writeptr, rdelta2);

printf("R1: ");
element_out_str(stdout, 0, R1);
printf("\n");
printf("R2: ");
element_out_str(stdout, 0, R2);
printf("\n");
printf("R3: ");
element_out_str(stdout, 0, R3);
printf("\n");
printf("R4: ");
element_out_str(stdout, 0, R4);
printf("\n");
printf("R5: ");
element_out_str(stdout, 0, R5);
printf("\n");

    element_clear(T1);
    element_clear(T2);
    element_clear(T3);
    element_clear(R1);
    element_clear(R2);
    element_clear(R3);
    element_clear(R4);
    element_clear(R5);
    element_clear(alpha);
    element_clear(beta);
    element_clear(c);
    element_clear(ralpha);
    element_clear(rbeta);
    element_clear(rx);
    element_clear(rdelta1);
    element_clear(rdelta2);
    //clear temp variables
    element_clear(z0);
    element_clear(z1);
    element_clear(e10);
    element_clear(et0);
}

int bbs_verify(unsigned char *sig,
	int hashlen, void *hash,
	bbs_group_public_key_t gpk)
{
    bbs_sys_param_ptr param = gpk->param;
    pairing_ptr pairing = param->pairing;
    field_ptr Fp = pairing->Zr;
    element_t T1, T2, T3;
    element_t R1, R2, R3, R4, R5;
    element_t c, salpha, sbeta, sx, sdelta1, sdelta2;
    element_t e10, e20, e21, et0, z0, z1;
    unsigned char *readptr = sig;
    UNUSED_VAR (hashlen);
    UNUSED_VAR (hash);

    element_init(T1, pairing->G1);
    element_init(T2, pairing->G1);
    element_init(T3, pairing->G1);
    element_init(R1, pairing->G1);
    element_init(R2, pairing->G1);
    element_init(R3, pairing->GT);
    element_init(R4, pairing->G1);
    element_init(R5, pairing->G1);

    element_init(c, Fp);
    element_init(salpha, Fp);
    element_init(sbeta, Fp);
    element_init(sx, Fp);
    element_init(sdelta1, Fp);
    element_init(sdelta2, Fp);

    element_init(e10, pairing->G1);
    element_init(e20, pairing->G2);
    element_init(e21, pairing->G2);
    element_init(et0, pairing->GT);
    element_init(z0, Fp);
    element_init(z1, Fp);

    readptr += element_from_bytes(T1, readptr);
    readptr += element_from_bytes(T2, readptr);
    readptr += element_from_bytes(T3, readptr);
    readptr += element_from_bytes(c, readptr);
    readptr += element_from_bytes(salpha, readptr);
    readptr += element_from_bytes(sbeta, readptr);
    readptr += element_from_bytes(sx, readptr);
    readptr += element_from_bytes(sdelta1, readptr);
    readptr += element_from_bytes(sdelta2, readptr);

    element_neg(z0, c);

    //element_pow(R1, gpk->u, to_mpz(salpha));
    //element_pow(e10, T1, to_mpz(z0));
    //element_mul(R1, R1, e10);
    element_pow2(R1, gpk->u, to_mpz(salpha), T1, to_mpz(z0));

    //element_pow(R2, gpk->v, to_mpz(sbeta));
    //element_pow(e10, T2, to_mpz(z0));
    //element_mul(R2, R2, e10);
    element_pow2(R2, gpk->v, to_mpz(sbeta), T2, to_mpz(z0));

    element_neg(z0, sdelta1);
    //element_pow(R4, gpk->u, to_mpz(z0));
    //element_pow(e10, T1, to_mpz(sx));
    //element_mul(R4, R4, e10);
    element_pow2(R4, gpk->u, to_mpz(z0), T1, to_mpz(sx));

    element_neg(z0, sdelta2);
    //element_pow(R5, gpk->v, to_mpz(z0));
    //element_pow(e10, T2, to_mpz(sx));
    //element_mul(R5, R5, e10);
    element_pow2(R5, gpk->v, to_mpz(z0), T2, to_mpz(sx));


    /*
     * compute R3 more efficiently.  use precomputed e(g1,g2)^{-1},
     * e(h,g2), and e(h,w).  this leaves e(T3,g2)^sx and e(T3,w)^c;
     * compute these with one pairing as e(T3, g2^sx w^c).
     */

    //element_pow(e20, gpk->g2, to_mpz(sx));
    //element_pow(e21, gpk->w, to_mpz(c));
    //element_mul(e20, e20, e21);
    element_pow2(e20, gpk->g2, to_mpz(sx), gpk->w, to_mpz(c));
    bilinear_map(R3, T3, e20, pairing);

    //element_pow(et0, gpk->pr_g1_g2_inv, to_mpz(c));
    //element_mul(R3, R3, et0);

    element_add(z0, salpha, sbeta);
    element_neg(z0, z0);
    //element_pow(et0, gpk->pr_h_w, to_mpz(z0));
    //element_mul(R3, R3, et0);

    element_add(z1, sdelta1, sdelta2);
    element_neg(z1, z1);
    //element_pow(et0, gpk->pr_h_g2, to_mpz(z1));

    element_pow3(et0, gpk->pr_g1_g2_inv, to_mpz(c),
                 gpk->pr_h_w, to_mpz(z0), gpk->pr_h_g2, to_mpz(z1));
    element_mul(R3, R3, et0);

printf("R1: ");
element_out_str(stdout, 0, R1);
printf("\n");
printf("R2: ");
element_out_str(stdout, 0, R2);
printf("\n");
printf("R3: ");
element_out_str(stdout, 0, R3);
printf("\n");
printf("R4: ");
element_out_str(stdout, 0, R4);
printf("\n");
printf("R5: ");
element_out_str(stdout, 0, R5);
printf("\n");

    element_clear(T1);
    element_clear(T2);
    element_clear(T3);
    element_clear(R1);
    element_clear(R2);
    element_clear(R3);
    element_clear(R4);
    element_clear(R5);
    element_clear(c);
    element_clear(salpha);
    element_clear(sbeta);
    element_clear(sx);
    element_clear(sdelta1);
    element_clear(sdelta2);
    element_clear(e10);
    element_clear(et0);
    element_clear(z0);
    element_clear(z1);
    return 1;
}

int bbs_open(element_t A, bbs_group_public_key_t gpk, bbs_manager_private_key_t gmsk,
	int hashlen, void *hash, unsigned char *sig)
{
    bbs_sys_param_ptr param = gpk->param;
    pairing_ptr pairing = param->pairing;
    field_ptr Fp = pairing->Zr;
    element_t T1, T2, T3;
    element_t R1, R2, R3, R4, R5;
    element_t c, salpha, sbeta, sx, sdelta1, sdelta2;
    element_t e10, et0, z0;
    unsigned char *readptr = sig;
    int result;
    UNUSED_VAR (hashlen);
    UNUSED_VAR (hash);

    //TODO: consolidate with verify
    element_init(T1, pairing->G1);
    element_init(T2, pairing->G1);
    element_init(T3, pairing->G1);
    element_init(R1, pairing->G1);
    element_init(R2, pairing->G1);
    element_init(R3, pairing->GT);
    element_init(R4, pairing->G1);
    element_init(R5, pairing->G1);

    element_init(c, Fp);
    element_init(salpha, Fp);
    element_init(sbeta, Fp);
    element_init(sx, Fp);
    element_init(sdelta1, Fp);
    element_init(sdelta2, Fp);

    element_init(e10, pairing->G1);
    element_init(et0, pairing->GT);
    element_init(z0, Fp);

    readptr += element_from_bytes(T1, readptr);
    readptr += element_from_bytes(T2, readptr);
    readptr += element_from_bytes(T3, readptr);
    readptr += element_from_bytes(c, readptr);
    readptr += element_from_bytes(salpha, readptr);
    readptr += element_from_bytes(sbeta, readptr);
    readptr += element_from_bytes(sx, readptr);
    readptr += element_from_bytes(sdelta1, readptr);
    readptr += element_from_bytes(sdelta2, readptr);

    element_neg(z0, c);
    element_pow(R1, gpk->u, to_mpz(salpha));
    element_pow(e10, T1, to_mpz(z0));
    element_mul(R1, R1, e10);

    element_pow(R2, gpk->v, to_mpz(sbeta));
    element_pow(e10, T2, to_mpz(z0));
    element_mul(R2, R2, e10);

    element_neg(z0, sdelta1);
    element_pow(R4, gpk->u, to_mpz(z0));
    element_pow(e10, T1, to_mpz(sx));
    element_mul(R4, R4, e10);

    element_neg(z0, sdelta2);
    element_pow(R5, gpk->v, to_mpz(z0));
    element_pow(e10, T2, to_mpz(sx));
    element_mul(R5, R5, e10);

    bilinear_map(R3, T3, gpk->w, pairing);
    bilinear_map(et0, gpk->g1, gpk->g2, pairing);
    element_invert(et0, et0);
    element_mul(R3, R3, et0);
    element_pow(R3, R3, to_mpz(c));

    bilinear_map(et0, T3, gpk->g2, pairing);
    element_pow(et0, et0, to_mpz(sx));
    element_mul(R3, R3, et0);

    element_add(z0, salpha, sbeta);
    element_neg(z0, z0);
    bilinear_map(et0, gpk->h, gpk->w, pairing);
    element_pow(et0, et0, to_mpz(z0));
    element_mul(R3, R3, et0);

    element_add(z0, sdelta1, sdelta2);
    element_neg(z0, z0);
    bilinear_map(et0, gpk->h, gpk->g2, pairing);
    element_pow(et0, et0, to_mpz(z0));
    element_mul(R3, R3, et0);

    //if mismatch result = 0;
    //} else {

    element_pow(A, T1, gmsk->xi1);
    element_pow(e10, T2, gmsk->xi2);
    element_mul(A, A, e10);
    element_invert(A, A);
    element_mul(A, A, T3);
    result =1;
    //}

    element_clear(T1);
    element_clear(T2);
    element_clear(T3);
    element_clear(R1);
    element_clear(R2);
    element_clear(R3);
    element_clear(R4);
    element_clear(R5);
    element_clear(c);
    element_clear(salpha);
    element_clear(sbeta);
    element_clear(sx);
    element_clear(sdelta1);
    element_clear(sdelta2);
    element_clear(e10);
    element_clear(et0);
    element_clear(z0);

    return result;
}

int main(void)
{
    bbs_sys_param_t sp;
    bbs_group_public_key_t gpk;
    bbs_manager_private_key_t gmsk;
    bbs_group_private_key_t gsk[5];
    pairing_t pairing;
    unsigned char *sig;
    int result;
    element_t A;
    double t0, t1;

    pairing_init_inp_str(pairing, stdin);

    printf("gen sys param...\n");
    bbs_gen_sys_param(sp, pairing);
    printf("gen keys...\n");
    t0 = get_time();
    bbs_gen(gpk, gmsk, 5, gsk, sp);
    t1 = get_time();
    printf("%fs elapsed\n", t1 - t0);
    t0 = t1;
    printf("sign...\n");
    sig = (unsigned char *) malloc(sp->signature_length);
    bbs_sign(sig, 0, NULL, gpk, gsk[0]);
    t1 = get_time();
    printf("%fs elapsed\n", t1 - t0);
    t0 = t1;
    printf("verify...\n");
    result = bbs_verify(sig, 0, NULL, gpk);
    if (result) {
	printf("signature verifies\n");
    } else {
	printf("signature does not verify\n");
    }
    t1 = get_time();
    printf("%fs elapsed\n", t1 - t0);
    t0 = t1;
    element_init(A, pairing->G1);
    bbs_open(A, gpk, gmsk, 0, NULL, sig);
    printf("open A = ");
    element_out_str(stdout, 0, A);
    printf("\n");
    printf("gsk0 A = ");
    element_out_str(stdout, 0, gsk[0]->A);
    printf("\n");
    t1 = get_time();
    printf("%fs elapsed\n", t1 - t0);
    t0 = t1;

    return 0;
}
