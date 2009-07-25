// Input/output utility routines for pairing parameters.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_symtab.h"
#include "pbc_memory.h"
#include "param_io.h"

void param_out_type(FILE *stream, char *s) {
  fprintf(stream, "type %s\n", s);
}

void param_out_mpz(FILE *stream, char *s, mpz_t z) {
  fprintf(stream, "%s ", s);
  mpz_out_str(stream, 0, z);
  fprintf(stream, "\n");
}

void param_out_int(FILE *stream, char *s, int i) {
  mpz_t z;
  mpz_init(z);

  mpz_set_si(z, i);
  param_out_mpz(stream, s, z);
  mpz_clear(z);
}

void lookup_mpz(mpz_t z, const char *(*tab)(const char *), const char *key) {
  const char *data = tab(key);
  if (!data) {
    pbc_error("missing param: `%s'", key);
    return;
  }
  mpz_set_str(z, data, 0);
}

int lookup_int(const char *(*tab)(const char *), const char *key) {
  int res;
  mpz_t z;
  const char *data = tab(key);
  if (!data) {
    pbc_error("missing param: `%s'", key);
    return 0;
  }

  mpz_init(z);

  mpz_set_str(z, data, 0);
  res = mpz_get_si(z);

  mpz_clear(z);
  return res;
}
