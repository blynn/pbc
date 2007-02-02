#include "pbc.h"

int main(int argc, char **argv)
{
    darray_t L;
    unsigned int D = 0;
    cm_info_ptr cm;
    darray_init(L);

    if (argc > 1) {
	D = atoi(argv[1]);
	if (D % 120) {
	    fprintf(stderr, "D must be multiple of 120\n");
	    exit(1);
	}
    }

    void try(int tryD) {
	int qbits, rbits;

	if (find_freeman_curve(L, tryD, 500)) {
	    while (darray_count(L)) {
		cm = darray_at(L, 0);
		qbits = mpz_sizeinbase(cm->q, 2);
		rbits = mpz_sizeinbase(cm->r, 2);
		printf("%d, %d, %d\n", tryD, qbits, rbits);
		fflush(stdout);
		darray_remove_index(L, 0);
		cm_info_clear(cm);
	    }
	}
    }

    printf("D < %u, bits in q, bits in r\n", 1000000000);
    while (D < 1000000000) {
	try(D + 43);
	try(D + 67);
	D+=120;
    }

    return 0;
}
