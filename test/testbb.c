//Boneh-Boyen signatures demo
#include "pbc_sig.h"

int main(void)
{
    pairing_t pairing;
    bb_sys_param_t param;
    bb_public_key_t pk;
    bb_private_key_t sk;
    unsigned char *sig;

    printf("initializing pairing from stdin...\n");

    pairing_init_inp_str(pairing, stdin);
    printf("generating BB system parameters...\n");
    bb_gen_sys_param(param, pairing);
    printf("generating key pair...\n");
    bb_gen(pk, sk, param);

    sig = (unsigned char *) malloc(param->signature_length);

    printf("signing...\n");
    bb_sign(sig, 11, (unsigned char *) "hello world", pk, sk);

    printf("verifying...\n");
    if (bb_verify(sig, 11, (unsigned char *) "hello world", pk)) {
	printf("signature verifies\n");
    } else {
	printf("signature does not verify\n");
    }
    return 0;
}
