#include "pbc.h"

int main(int argc, char **argv)
{
    int rbits = 160, qbits = 512;
    e_param_t ep;

    e_param_init(ep);

    if (argc > 1) {
	rbits = atoi(argv[1]);
    }
    if (argc > 2) {
	qbits = atoi(argv[2]);
    }

    e_param_gen(ep, rbits, qbits);
    e_param_out_str(stdout, ep);

    return 0;
}
