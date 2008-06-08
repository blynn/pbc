// Generate MNT curve(s) for a given D.

#include "pbc.h"
#include "pbc_utils.h"  // for UNUSED_VAR

void generate(cm_info_t cm, void *data)
{
    UNUSED_VAR(data);
    d_param_t param;
    d_param_init(param);

    fprintf(stderr, "gendparam: computing Hilbert polynomial and finding roots...\n");
    d_param_from_cm(param, cm);
    fprintf(stderr, "gendparam: bits in q = %zu\n", mpz_sizeinbase(cm->q, 2));
    d_param_out_str(stdout, param);
}

int main(int argc, char **argv)
{
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

    if (! find_mnt6_curve(generate, NULL, D, 500)) {
	fprintf(stderr, "No suitable curves for this D\n");
	exit(1);
    }
    return 0;
}
