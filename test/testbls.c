//Boneh-Lynn-Shacham short signatures demo
//This one uses the signature library and is more realistic
#include "pbc_sig.h"

int main(void)
{
    pairing_t pairing;
    bls_sys_param_t param;
    bls_public_key_t pk;
    bls_private_key_t sk;
    unsigned char *sig;

    printf("reading pairing from stdin...\n");

    pairing_init_inp_str(pairing, stdin);
    printf("generating BLS system parameters...\n");
    bls_gen_sys_param(param, pairing);
    printf("generating key pair...\n");
    bls_gen(pk, sk, param);

    sig = (unsigned char *) malloc(param->signature_length);

    printf("signing...\n");
    bls_sign(sig, 11, (unsigned char *) "hello world", sk);

    printf("verifying...\n");
    if (bls_verify(sig, 11, (unsigned char *) "hello world", pk)) {
	printf("signature verifies\n");
    } else {
	printf("signature does not verify\n");
    }
    return 0;
}
