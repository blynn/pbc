#include "pbc.h"

int main(void)
{
    element_t g, h, x;
    element_t zg, zh;
    pairing_t pairing;
    mpz_t z;

    pairing_init_inp_str(pairing, stdin);

    element_init_G1(g, pairing);
    element_init_G1(zg, pairing);
    element_init_G2(h, pairing);
    element_init_G2(zh, pairing);
    element_init_GT(x, pairing);
    element_random(g);
    element_random(h);
    //pairing->phi(g, h, pairing);
    element_printf("g = %B\n", g);
    element_printf("h = %B\n", h);
    bilinear_map(x, g, h, pairing);
    element_printf("f(g, h) = %B\n", x);

    mpz_init(z);
    pbc_mpz_random(z, pairing->r);
    element_printf("z = %Z\n", z);

    element_pow(x, x, z);
    element_printf("f(g, h)^z = %B\n", x);

    element_pow(zg, g, z);
    element_printf("g^z = %B\n", zg);
    bilinear_map(x, zg, h, pairing);
    element_printf("f(g^z, h) = %B\n", x);

    element_pow(zh, h, z);
    element_printf("h^z = %B\n", zh);
    bilinear_map(x, g, zh, pairing);
    element_printf("f(g, h^z) = %B\n", x);

    {
	int i;
	int len = element_length_in_bytes(h);
	printf("length_in_bytes(h) = %d\n", len);
	unsigned char *data = malloc(len);
	element_to_bytes(data, h);
	for (i=0; i<len; i++) {
	    printf(" %02X", data[i]);
	    if (15 == (i % 16)) printf("\n");
	}
	printf("\n");
	element_from_bytes(h, data);
	element_printf("from_bytes h = %B\n", h);
    }
    return 0;
}
