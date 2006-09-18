#include "pbc.h"

int main(void)
{
    element_t g, h, x;
    element_t zg, zh;
    pairing_t pairing;
    mpz_t z;

    pairing_init_inp_str(pairing, stdin);

    element_init(g, pairing->G1);
    element_init(zg, pairing->G1);
    element_init(h, pairing->G2);
    element_init(zh, pairing->G2);
    element_init(x, pairing->GT);
    element_random(g);
    element_random(h);
    //pairing->phi(g, h, pairing);
    printf("g = ");
    element_out_str(stdout, 0, g);
    printf("\n");
    printf("h = ");
    element_out_str(stdout, 0, h);
    printf("\n");
    bilinear_map(x, g, h, pairing);
    printf("f(g, h) = ");
    element_out_str(stdout, 0, x);
    printf("\n");

    mpz_init(z);
    pbc_mpz_random(z, pairing->r);
    printf("z = ");
    mpz_out_str(stdout, 0, z);
    printf("\n");

    element_pow(x, x, z);
    printf("f(g, h)^z = ");
    element_out_str(stdout, 0, x);
    printf("\n");

    element_pow(zg, g, z);
    printf("g^z = ");
    element_out_str(stdout, 0, zg);
    printf("\n");
    bilinear_map(x, zg, h, pairing);
    printf("f(g^z, h) = ");
    element_out_str(stdout, 0, x);
    printf("\n");

    element_pow(zh, h, z);
    printf("h^z = ");
    element_out_str(stdout, 0, zh);
    printf("\n");
    bilinear_map(x, g, zh, pairing);
    printf("f(g, h^z) = ");
    element_out_str(stdout, 0, x);
    printf("\n");

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
	printf("from_bytes h = ");
	element_out_str(stdout, 0, h);
	printf("\n");
    }
    return 0;
}

/*
int main(void)
{
    element_t x, x2, y, y2, r;
    pairing_t pairing;

    //gmp_leak_check();
    pairing_init_inp_str(pairing, stdin);

    element_init(x, pairing->G1);
    element_init(y, pairing->G2);
    element_init(x2, pairing->G1);
    element_init(y2, pairing->G2);
    element_init(r, pairing->GT);

    if (1) {
	int i;
	mpz_t pow;
	mpz_init(pow);
	element_random(y);
	element_random(y2);
	for (i=0; i<10; i++) {
	pbc_mpz_random(pow, pairing->r);
	    element_pow(y2, y, pow);
	    //element_pow(x2, x, pow);
	}
	mpz_clear(pow);
	mem_report();
    } else {
	element_random(x);
	element_random(y);
	printf("x = ");
	element_out_str(stdout, 0, x);
	printf("\n");
	printf("y = ");
	element_out_str(stdout, 0, y);
	printf("\n");
	element_mul(x2, x, x);
	element_mul(y2, y, y);

	bilinear_map(r, x2, y, pairing);
	printf("e(x^2,y) = ");
	element_out_str(stdout, 0, r);
	printf("\n");

	bilinear_map(r, x, y2, pairing);
	printf("e(x,y^2) = ");
	element_out_str(stdout, 0, r);
	printf("\n");

	bilinear_map(r, x, y, pairing);
	element_mul(r, r, r);
	printf("e(x,y)^2 = ");
	element_out_str(stdout, 0, r);
	printf("\n");
    }

    return 0;
}
*/
