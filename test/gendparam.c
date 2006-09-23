#include "pbc.h"
#include "d_param.h"

int main(int argc, char **argv)
{
    darray_t L;
    d_param_t param;
    cm_info_ptr cm;
    int D = 9563;

    if (argc > 1) {
	int m;
	D = atoi(argv[1]);
	m = D % 4;
	if (D <= 0 || m == 1 || m == 2) {
	    fprintf(stderr, "D must be 0 or 3 mod 4 and positive\n");
	    exit(1);
	}
    }
    fprintf(stderr, "Using D = %d\n", D);
    darray_init(L);

    find_mnt6_curve(L, D, 500);

    if (!L->count) {
	fprintf(stderr, "No suitable curves for this D\n");
	exit(1);
    }
    cm = darray_at(L, 0);
    d_param_init(param);

    fprintf(stderr, "gencparam: computing Hilbert polynomial and finding roots...\n");
    d_param_from_cm(param, cm);
    fprintf(stderr, "gencparam: bits in q = %zu\n", mpz_sizeinbase(cm->q, 2));
    d_param_out_str(stdout, param);
    return 0;
}
