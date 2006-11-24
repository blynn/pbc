#include "pbc.h"
#include "pbc_fp.h"

int main(void)
{
    field_t zp, rx, fp2;
    mpz_t prime;
    element_t a;
    element_t f, g;

    mpz_init(prime);
    mpz_set_ui(prime, 100000);
    mpz_setbit(prime, 33);
    mpz_nextprime(prime, prime);
    element_printf("prime = %Z\n", prime);

    field_init_fp(zp, prime);
    field_init_poly(rx, zp);
    element_init(f, rx);
    element_init(g, rx);
    element_init(a, zp);

    element_set_si(a, 1);
    poly_set_coeff(f, a, 2);
    poly_set_coeff(f, a, 0);

    element_out_str(stdout, 0, f);
    if (poly_is_irred(f)) {
	printf(" is irreducible\n");
    } else {
	printf(" is not irreducible\n");
    }

    {
	unsigned char *data;
	int i;
	int n = element_length_in_bytes(f);
	printf("serialized f =");
	data = (unsigned char *) pbc_malloc(n);
	element_to_bytes(data, f);
	for (i=0; i<n; i++) {
	    printf(" %02X", data[i]);
	}
	printf("\n");
	printf("deserialize check = ");
	element_from_bytes(f, data);
	element_out_str(stdout, 10, f);
	printf("\n");
    }

    field_init_polymod(fp2, f);
    element_clear(f);
    element_clear(g);
    element_init(f, fp2);
    element_init(g, fp2);
    element_random(f);
    element_random(g);
    element_printf("f: %B, g: %B\n", f, g);
    element_invert(g, f);
    element_printf("inv f: %B\n", g);
    element_mul(g, f, g);
    element_printf("prod: %B\n", g);
    {
	do {
	    element_random(f);
	} while (!element_is_sqr(f));
	element_printf("random square f: %B\n", f);
	element_sqrt(f, f);
	element_printf("sqrt f: %B\n", f);
	element_mul(f, f, f);
	element_printf("f: %B\n", f);
    }
    {
	unsigned char *data;
	int i;
	int n = element_length_in_bytes(f);
	printf("serialized f =");
	data = (unsigned char *) pbc_malloc(n);
	element_to_bytes(data, f);
	for (i=0; i<n; i++) {
	    printf(" %02X", data[i]);
	}
	printf("\n");
	printf("deserialize check = ");
	element_random(f);
	element_from_bytes(f, data);
	element_out_str(stdout, 10, f);
	printf("\n");
    }

    element_clear(f);
    element_clear(g);
    return 0;
}
