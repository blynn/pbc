#include "pbc.h"

int main(int argc, char **argv)
{
    darray_t L;
    unsigned int D = 7;
    cm_info_ptr cm;
    int i=0;
    darray_init(L);

    if (argc > 1) {
	D = atoi(argv[1]);
	if (D < 7 || (D % 4) != 3) {
	    fprintf(stderr, "D must be 3 mod 4 and at least 7\n");
	    exit(1);
	}
    }

    void try(void) {
	int qbits, rbits;

	if (find_mnt6_curve(L, D, 500)) {
	    for (; i<L->count; i++) {
		cm = L->item[i];
		qbits = mpz_sizeinbase(cm->q, 2);
		rbits = mpz_sizeinbase(cm->r, 2);
		printf("%d, %d, %d\n", D, qbits, rbits);
		fflush(stdout);
	    }
	}
    }

    printf("D < %u, bits in q, bits in r\n", 1000000000);
    while (D < 1000000000) {
	try();
	D++;
	try();
	D+=3;
    }

    return 0;
}
