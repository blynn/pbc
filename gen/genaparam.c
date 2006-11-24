#include <stdlib.h>
#include "pbc.h"

int main(int argc, char **argv)
{
    int rbits = 160, qbits = 512;
    a_param_t ap;

    a_param_init(ap);

    if (argc > 1) {
	rbits = atoi(argv[1]);
    }
    if (argc > 2) {
	qbits = atoi(argv[2]);
    }
    a_param_gen(ap, rbits, qbits);

    a_param_out_str(stdout, ap);
    a_param_clear(ap);

    return 0;
}
