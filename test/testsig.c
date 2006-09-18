//Boneh-Lynn-Shacham short signatures demo
//see sig.c, testbls.c for a more practical implementation
#include "pbc.h"

int main(void)
{
    element_t g, h, s;
    element_t pubkey, sig;
    pairing_t pairing;
    mpz_t secret;

    pairing_init_inp_str(pairing, stdin);

    element_init(g, pairing->G2);
    element_init(pubkey, pairing->G2);
    element_init(h, pairing->G1);
    element_init(sig, pairing->G1);
    element_init(s, pairing->GT);
    mpz_init(secret);

    printf("Short signature test\n");

    //generate system parameters
    element_random(g);
    printf("system parameter g = ");
    element_out_str(stdout, 0, g);
    printf("\n");

    //generate private key
    pbc_mpz_random(secret, pairing->r);
    printf("private key = ");
    mpz_out_str(stdout, 0, secret);
    printf("\n");

    //compute corresponding public key
    element_pow(pubkey, g, secret);
    printf("public key = ");
    element_out_str(stdout, 0, pubkey);
    printf("\n");

    //generate element from a hash
    //for toy examples, should check that pairing(g, h) != 1
    element_from_hash(h, 13, "hashofmessage");
    printf("message hash = ");
    element_out_str(stdout, 0, h);
    printf("\n");

    //h^secret is the signature
    //in real life: only output the first coordinate
    element_pow(sig, h, secret);
    printf("signature = ");
    element_out_str(stdout, 0, sig);
    printf("\n");

    {
	int n = element_length_in_bytes_compressed(sig);
	int i;
	unsigned char *data = malloc(n);
	element_to_bytes_compressed(data, sig);
	printf("compressed = ");
	for (i=0; i<n; i++) {
	    printf("%02X", data[i]);
	}
	printf("\n");

	element_from_bytes_compressed(sig, data);
	printf("uncompressed = ");
	element_out_str(stdout, 0, sig);
	printf("\n");
    }

    //verification part 1
    bilinear_map(s, sig, g, pairing);
    printf("f(sig, g) = ");
    element_out_str(stdout, 0, s);
    printf("\n");

    //verification part 2
    //should match above
    bilinear_map(s, h, pubkey, pairing);
    printf("f(message hash, pubkey) = ");
    element_out_str(stdout, 0, s);
    printf("\n");

    return 0;
}
