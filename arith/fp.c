// F_p initialization.
//
// Specific implementations of F_p are found in naivefp.c, fastfp.c, fasterfp.c
// and montfp.c. For pairing-based cryptosystems, montfp.c is the fastest.
// I keep all versions around for testing, and also to show off the modularity
// of the code.

#include <stdarg.h>
#include <stdio.h>
#include <stdint.h> // for intptr_t
#include <gmp.h>
#include <string.h>
#include "pbc_utils.h"
#include "pbc_field.h"
#include "pbc_fp.h"

// By default, use the montfp.c implementation of F_p. After
// pbc_tweak_use_fp(), future field_init_fp calls will use the specified
// implementation. This is useful for benchmarking and testing.
static void (*option_fpinit) (field_ptr f, mpz_t prime) = field_init_mont_fp;

void pbc_tweak_use_fp(char *s) {
  if (!strcmp(s, "naive")) {
    option_fpinit = field_init_naive_fp;
  } else if (!strcmp(s, "fast")) {
    option_fpinit = field_init_fast_fp;
  } else if (!strcmp(s, "faster")) {
    option_fpinit = field_init_faster_fp;
  } else if (!strcmp(s, "mont")) {
    option_fpinit = field_init_mont_fp;
  } else {
    pbc_error("no such Fp implementation: %s", s);
  }
}

void field_init_fp(field_ptr f, mpz_t modulus) {
  if (mpz_fits_ulong_p(modulus)) {
    // If this case mattered, I'd have written a F_p implementation specialized
    // for moduli that fits into machine words.
    field_init_naive_fp(f, modulus);
  } else {
    if (mpz_odd_p(modulus)) {
      option_fpinit(f, modulus);
    } else {
      // montfp.c only supports odd moduli.
      field_init_faster_fp(f, modulus);
    }
  }
}
