#include "pbc.h"

int main(void)
{
    field_t fp, fp2;
    mpz_t prime;
    element_t a, b, c;

    mpz_init(prime);
    mpz_set_ui(prime, 82);
    mpz_nextprime(prime, prime);
    field_init_fp(fp, prime);
    field_init_fi(fp2, fp);
    element_init(a, fp2);
    element_init(b, fp2);
    element_init(c, fp2);

    printf("field: ");
    mpz_out_str(stdout, 0, prime);
    printf("^2\n");

    element_random(a);
    element_random(b);
    printf("a = ");
    element_out_str(stdout, 0, a);
    printf(", b = ");
    element_out_str(stdout, 0, b);
    printf("\n");

    element_add(c, a, b);
    printf("a + b = ");
    element_out_str(stdout, 0, c);
    printf("\n");

    element_mul(c, a, b);
    printf("a * b = ");
    element_out_str(stdout, 0, c);
    printf("\n");

    for (;;) {
	element_random(a);
	printf("new a = ");
	element_out_str(stdout, 0, a);
	printf("\n");

	if (element_is_sqr(a)) break;
	printf(" is not a square\n");
    }
    element_sqrt(c, a);
    printf("sqrt(a) = ");
    element_out_str(stdout, 0, c);
    printf("\n");
    element_mul(c, c, c);
    printf("sqrt(a) * sqrt(a) = ");
    element_out_str(stdout, 0, c);
    printf("\n");
    element_invert(c, a);
    printf("1/a = ");
    element_out_str(stdout, 0, c);
    printf("\n");
    element_mul(c, c, a);
    printf("1/a * a = ");
    element_out_str(stdout, 0, c);
    printf("\n");

    return 0;
}
