#ifndef RANDOM_H
#define RANDOM_H

#include <gmp.h>

void random_set_file(char *filename);
void random_set_deterministic(void);
void pbc_mpz_random(mpz_t z, mpz_t limit);
void random_push(void (*random_fn)(mpz_t, mpz_t, void *), void *random_data);
void random_pop(void);

#endif //RANDOM_H
