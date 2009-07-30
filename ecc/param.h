// Input/output routines common to all pairing parameters.

// Requires:
// * param.h
// * stdio.h
// * gmp.h
#ifndef __PARAM_UTILS_H__
#define __PARAM_UTILS_H__

#pragma GCC visibility push(hidden)

void param_out_type(FILE *stream, char *s);
void param_out_mpz(FILE *stream, char *s, mpz_t z);
void param_out_int(FILE *stream, char *s, int i);
// TODO: Replace with a stdarg function, e.g.
//   err = lookup("ZZi", "p", "n", "l",  p->p, p->n, &p->l);
int lookup_int(int *n, const char *(*tab)(const char *), const char *key);
int lookup_mpz(mpz_t z, const char *(*tab)(const char *), const char *key);

#pragma GCC visibility pop

#endif //__PARAM_UTILS_H__
