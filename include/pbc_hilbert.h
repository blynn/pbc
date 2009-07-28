// Requires:
// * gmp.h
#ifndef __PBC_HILBERT_H__
#define __PBC_HILBERT_H__

// Allocate an array of mpz_t and fill it with the coefficients of the Hilbert
// polynomial H_D(x). Returns the size of array.
size_t pbc_hilbert(mpz_t **arr, int D);

// Free an array allocated by `pbc_hilbert()`.
void pbc_hilbert_free(mpz_t *arr, size_t n);

#endif //__PBC_HILBERT_H__
