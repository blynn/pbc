//Boneh-Lynn-Shacham short signatures demo
//see sig.c, testbls.c for a more practical implementation
#include "pbc.h"

int main(void)
{
    element_t g, h, s;
    element_t pubkey, sig;
    pairing_t pairing;
    element_t secret;

    pairing_init_inp_str(pairing, stdin);

    element_init_G2(g, pairing);
    element_init_G2(pubkey, pairing);
    element_init_G1(h, pairing);
    element_init_G1(sig, pairing);
    element_init_GT(s, pairing);
    element_init_Zr(secret, pairing);

    printf("Short signature test\n");

    //generate system parameters
    element_random(g);
    element_printf("system parameter g = %B\n", g);

    //generate private key
    element_random(secret);
    element_printf("private key = %B\n", secret);

    //compute corresponding public key
    element_pow_zn(pubkey, g, secret);
    element_printf("public key = %B\n", pubkey);

    //generate element from a hash
    //for toy examples, should check that pairing(g, h) != 1
    element_from_hash(h, 13, "hashofmessage");
    element_printf("message hash = %B\n", h);

    //h^secret is the signature
    //in real life: only output the first coordinate
    element_pow_zn(sig, h, secret);
    element_printf("signature = %B\n", sig);

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
	element_printf("uncompressed = %B\n", sig);
    }

    //verification part 1
    bilinear_map(s, sig, g, pairing);
    element_printf("f(sig, g) = %B\n", s);

    //verification part 2
    //should match above
    bilinear_map(s, h, pubkey, pairing);
    element_printf("f(message hash, pubkey) = %B\n", s);

    return 0;
}
