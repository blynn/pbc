#include "pbc.h"

int main(void)
{
    element_t g, h, x;
    element_t zg, zh, z;
    pairing_t pairing;

    pairing_init_inp_str(pairing, stdin);

    element_init_G1(g, pairing);
    element_init_G1(zg, pairing);
    element_init_G2(h, pairing);
    element_init_G2(zh, pairing);
    element_init_GT(x, pairing);
    element_init_Zr(z, pairing);
    element_random(g);
    element_random(h);
    //pairing->phi(g, h, pairing);
    element_printf("g = %B\n", g);
    element_printf("h = %B\n", h);
    bilinear_map(x, g, h, pairing);
    element_printf("f(g, h) = %B\n", x);

    element_random(z);
    element_printf("z = %B\n", z);

    element_pow_zn(x, x, z);
    element_printf("f(g, h)^z = %B\n", x);

    element_pow_zn(zg, g, z);
    element_printf("g^z = %B\n", zg);
    bilinear_map(x, zg, h, pairing);
    element_printf("f(g^z, h) = %B\n", x);

    element_pow_zn(zh, h, z);
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
