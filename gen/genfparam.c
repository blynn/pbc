#include "pbc.h"

int main(int argc, char **argv)
{
    int bits = 160;
    if (argc > 1) {
	bits = atoi(argv[1]);
	if (bits < 1) {
	    fprintf(stderr, "Usage: %s [BITS]\n", argv[0]);
	    exit(1);
	}
    }
    f_param_t fp;

    f_param_init(fp);
    f_param_gen(fp, bits);
    f_param_out_str(stdout, fp);

    return 0;
}
