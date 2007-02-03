#include "pbc.h"

int main(int argc, char **argv)
{
    darray_t L;
    g_param_t param;
    cm_info_ptr cm;
    int D = 35707;

    if (argc > 1) {
	int m;
	D = atoi(argv[1]);
	m = D % 120;
	if (D <= 0 || (m != 43 && m != 67)) {
	    fprintf(stderr, "D must be 43 or 67 mod 120 and positive\n");
	    exit(1);
	}
    }
    fprintf(stderr, "Using D = %d\n", D);
    darray_init(L);

    find_freeman_curve(L, D, 500);

    if (!L->count) {
	fprintf(stderr, "No suitable curves for this D\n");
	exit(1);
    }
    cm = darray_at(L, 0);
    g_param_init(param);

    fprintf(stderr, "gengparam: computing Hilbert polynomial and finding roots...\n");
    g_param_from_cm(param, cm);
    fprintf(stderr, "gengparam: bits in q = %zu\n", mpz_sizeinbase(cm->q, 2));
    g_param_out_str(stdout, param);

    // If we weren't exiting now, would cm_info_clear every entry of L
    // and darray_clear L
    return 0;
}
