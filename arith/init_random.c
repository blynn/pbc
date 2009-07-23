#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "pbc_assert.h"
#include "pbc_random.h"

void init_random_function(void) {
  FILE *fp;
  fp = fopen("/dev/urandom", "rb");
  if (!fp) {
    pbc_warn("could not open /dev/urandom, using deterministic random number generator");
    random_set_deterministic();
  } else {
    random_set_file("/dev/urandom");
    fclose(fp);
  }
}
