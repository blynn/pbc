//Cha-Cheon and Sakai-Kasahara-Schnorr Identity-Based Signatures demo
#include "sig.h"

int main(void)
{
    pairing_t pairing;
    ib_sys_param_t param;
    ib_master_key_t master;
    ib_private_key_t sk;
    unsigned char *sig;

    printf("reading pairing from stdin...\n");

    pairing_init_inp_str(pairing, stdin);
    printf("generating identity-based system parameters...\n");
    ib_setup(param, master, pairing);
    printf("extracting private key...\n");
    ib_extract(sk, 5, (unsigned char *) "alice", master);

    sig = (unsigned char *) malloc(cc_signature_length(param));

    printf("Cha-Cheon signatures:\n");
    printf("signing...\n");
    cc_sign(sig, 11, (unsigned char *) "hello world", sk);

    {
	int i;
	int n = cc_signature_length(param);
	printf("signature: ");
	for (i=0; i<n; i++) {
	    printf("%02X", (unsigned int) sig[i]);
	}
	printf("\n");
    }

    printf("verifying...\n");
    if (cc_verify(sig, 11, (unsigned char *) "hello world",
		5, (unsigned char *) "alice", param)) {
	printf("signature verifies\n");
    } else {
	printf("signature does not verify\n");
    }

    free(sig);

    /* Cannot use this scheme: it is patented
    sig = (unsigned char *) malloc(skschnorr_signature_length(param));

    printf("Sakai-Kasahara-Schnorr signatures:\n");
    printf("signing...\n");
    skschnorr_sign(sig, 11, (unsigned char *) "hello world", sk);

    {
	int i;
	int n = skschnorr_signature_length(param);
	printf("signature: ");
	for (i=0; i<n; i++) {
	    printf("%02X", (unsigned int) sig[i]);
	}
	printf("\n");
    }

    printf("verifying...\n");
    if (skschnorr_verify(sig, 11, (unsigned char *) "hello world",
		5, (unsigned char *) "alice", param)) {
	printf("signature verifies\n");
    } else {
	printf("signature does not verify\n");
    }
    free(sig);
    */

    return 0;
}
