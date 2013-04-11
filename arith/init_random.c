#include <stdio.h>
#include <stdint.h> // for intptr_t
#include <stdlib.h>
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_random.h"

void pbc_init_random(void) {
  FILE *fp;
  fp = fopen("/dev/urandom", "rb");
  if (!fp) {
    pbc_warn("could not open /dev/urandom, using deterministic random number generator");
    pbc_random_set_deterministic(0);
  } else {
    pbc_random_set_file("/dev/urandom");
    fclose(fp);
  }
}
