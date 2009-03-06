#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "pbc_random.h"

void init_random_function() {
    FILE *fp;
    fp = fopen("/dev/urandom", "rb");
    if (!fp) {
	fprintf(stderr, "Warning: could not open /dev/urandom, using deterministic random number generator\n");
	random_set_deterministic();
    } else {
	random_set_file("/dev/urandom");
	fclose(fp);
    }
}
