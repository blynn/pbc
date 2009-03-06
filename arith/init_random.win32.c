// Win32 Compatibility Code added by Yulian Kalev.
#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <wincrypt.h>
#include <gmp.h>
#include "pbc_random.h"
#include "pbc_utils.h"
#include "pbc_memory.h"

static void win32_mpz_random(mpz_t r, mpz_t limit, void *data)
{
    UNUSED_VAR (data);
    HCRYPTPROV phProv;
    if (!CryptAcquireContext(&phProv,NULL,NULL,PROV_RSA_FULL,0))
    {
	fprintf(stderr,"Couldn't create CryptContext: %x\n",(int)GetLastError());
	return;
    }
    int n, bytecount, leftover;
    unsigned char *bytes;
    mpz_t z;
    mpz_init(z);
    n = mpz_sizeinbase(limit, 2);
    bytecount = (n + 7) / 8;
    leftover = n % 8;
    bytes = (unsigned char *) pbc_malloc(bytecount);
    for (;;) {
	CryptGenRandom(phProv,bytecount,(byte *)bytes);
	if (leftover) {
	    *bytes = *bytes % (1 << leftover);
	}
	mpz_import(z, bytecount, 1, 1, 0, 0, bytes);
	if (mpz_cmp(z, limit) < 0) break;
    }
    CryptReleaseContext(phProv,0);
    mpz_set(r, z);
    mpz_clear(z);
    pbc_free(bytes);
}

void init_random_function() {
    set_random_function(win32_mpz_random, NULL);
}
