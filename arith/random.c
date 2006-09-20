#include <stdio.h>
#include <stdlib.h>
#include "random.h"
#include "utils.h"

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
//TODO: generate some kind of warning on error
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
    bytes = (unsigned char *) malloc(bytecount);
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
}

static void (*current_mpz_random)(mpz_t, mpz_t, void *) = deterministic_mpz_random;
static void *current_random_data;

void pbc_mpz_random(mpz_t z, mpz_t limit)
{
    current_mpz_random(z, limit, current_random_data);
}

void random_set_deterministic(void)
{
    current_mpz_random = deterministic_mpz_random;
    current_random_data = NULL;
}

void random_set_file(char *filename)
{
    current_mpz_random = file_mpz_random;
    current_random_data = filename;
}

//TODO: use a linked list, or have some other way of managing RNGs

static void (*kludge_mpz_random)(mpz_t, mpz_t, void *) = deterministic_mpz_random;
static void *kludge_random_data;

void random_push(void (*random_fn)(mpz_t, mpz_t, void *), void *random_data)
{
    kludge_mpz_random = current_mpz_random;
    kludge_random_data = current_random_data;
    current_mpz_random = random_fn;
    current_random_data = random_data;
}

void random_pop()
{
    current_mpz_random = kludge_mpz_random;
    current_random_data = kludge_random_data;
}
