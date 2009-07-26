// Input/output routines common to all pairing parameters.

// Requires:
// * param.h
// * stdio.h
// * gmp.h
#ifndef __PARAM_UTILS_H__
#define __PARAM_UTILS_H__
void param_out_type(FILE *stream, char *s);
void param_out_mpz(FILE *stream, char *s, mpz_t z);
void param_out_int(FILE *stream, char *s, int i);
int lookup_int(const char *(*tab)(const char *), const char *key);
void lookup_mpz(mpz_t z, const char *(*tab)(const char *), const char *key);
#endif //__PARAM_UTILS_H__
