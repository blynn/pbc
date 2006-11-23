// requires
// * stdio.h
// * gmp.h
// * fops.h
// * symtab.h
#ifndef __PBC_PARAM_H__
#define __PBC_PARAM_H__

void param_out_type(FILE *stream, char *s);
void param_out_mpz(FILE *stream, char *s, mpz_t z);
void param_out_int(FILE *stream, char *s, int i);
void param_read_generic(symtab_t tab, fetch_ops_t fops, void *ctx);
void param_read_buf(symtab_t tab, const char *buf, size_t len);
void param_read_str(symtab_t tab, FILE *stream);
void param_clear_tab(symtab_t tab);
int lookup_int(symtab_t tab, char *key);
void lookup_mpz(mpz_t z, symtab_t tab, char *key);

#endif //__PBC_PARAM_H__
