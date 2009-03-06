#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "pbc_random.h"
#include "pbc_utils.h"
#include "pbc_memory.h"

static void deterministic_mpz_random(mpz_t z, mpz_t limit, void *data)
{

    static gmp_randstate_t rs;
    static int rs_is_ready;
    UNUSED_VAR (data);

    if (!rs_is_ready) {
	gmp_randinit_default(rs);
	rs_is_ready = 1;
    }
    mpz_urandomm(z, rs, limit);
}

static void file_mpz_random(mpz_t r, mpz_t limit, void *data)
// TODO: Warn on error.
{
    char *filename = (char *) data;
    FILE *fp;
    int n, bytecount, leftover;
    unsigned char *bytes;
    mpz_t z;
    mpz_init(z);
    fp = fopen(filename, "rb");
    if (!fp) return;
    n = mpz_sizeinbase(limit, 2);
    bytecount = (n + 7) / 8;
    leftover = n % 8;
    bytes = (unsigned char *) pbc_malloc(bytecount);
    for (;;) {
	fread(bytes, 1, bytecount, fp);
	if (leftover) {
	    *bytes = *bytes % (1 << leftover);
	}
	mpz_import(z, bytecount, 1, 1, 0, 0, bytes);
	if (mpz_cmp(z, limit) < 0) break;
    }
    fclose(fp);
    mpz_set(r, z);
    mpz_clear(z);
    pbc_free(bytes);
}

static void (*current_mpz_random)(mpz_t, mpz_t, void *);
static void *current_random_data;
static int random_function_ready = 0;

void set_random_function(void (*fun)(mpz_t, mpz_t, void *), void *data) {
    current_mpz_random = fun;
    current_random_data = data;
    random_function_ready = 1;
}

void pbc_mpz_random(mpz_t z, mpz_t limit)
{
    if (!random_function_ready) {
        init_random_function();
    }
    current_mpz_random(z, limit, current_random_data);
}

void pbc_mpz_randomb(mpz_t z, unsigned int bits)
{
    mpz_t limit;
    mpz_init(limit);
    mpz_setbit(limit, bits);
    pbc_mpz_random(z, limit);
    mpz_clear(limit);
}

void random_set_deterministic(void)
{
    set_random_function(deterministic_mpz_random, NULL);
}

void random_set_file(char *filename)
{
    set_random_function(file_mpz_random, filename);
}
