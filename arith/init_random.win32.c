// Win32 Compatibility Code added by Yulian Kalev and Stefan Georg Weber.
#include <stdio.h>
#include <stdint.h> // for intptr_t
#include <stdlib.h>
#include <windows.h>
#include <wincrypt.h>
#include <gmp.h>
#include "pbc_random.h"
#include "pbc_utils.h"
#include "pbc_memory.h"

static void win32_mpz_random(mpz_t r, mpz_t limit, void *data) {
  UNUSED_VAR (data);
  HCRYPTPROV phProv;
  unsigned int error;
  if (!CryptAcquireContext(&phProv,NULL,NULL,PROV_RSA_FULL,0)) {
    error = GetLastError();
    if (error == 0x80090016) { //need to create a new keyset
      if (!CryptAcquireContext(&phProv,NULL,NULL,PROV_RSA_FULL,CRYPT_NEWKEYSET)) {
        pbc_error("Couldn't create CryptContext: %x", (int)GetLastError());
        return;
      }
    } else {
      pbc_error("Couldn't create CryptContext: %x", error);
      return;
    }
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

void pbc_init_random(void) {
  pbc_random_set_function(win32_mpz_random, NULL);
}
