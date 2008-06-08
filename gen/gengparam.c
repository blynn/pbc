// Generate Freeman curves with a given discriminant.
#include "pbc.h"

void generate(cm_info_t cm, void *data)
{
    g_param_t param;
    g_param_init(param);

    fprintf(stderr, "gengparam: computing Hilbert polynomial and finding roots...\n");
    g_param_from_cm(param, cm);
    fprintf(stderr, "gengparam: bits in q = %zu\n", mpz_sizeinbase(cm->q, 2));
    g_param_out_str(stdout, param);
}

int main(int argc, char **argv)
{
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

    if (!find_freeman_curve(generate, NULL, D, 500)) {
	fprintf(stderr, "No suitable curves for this D\n");
	exit(1);
    }

    return 0;
}
