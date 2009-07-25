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

struct pairing_s;
struct pbc_param_interface_s {
  void (*clear)(void *);
  void (*init_pairing)(struct pairing_s *, void *);
  void (*out_str)(FILE *stream, void *data);
};
typedef struct pbc_param_interface_s pbc_param_interface_t[1];
typedef struct pbc_param_interface_s *pbc_param_interface_ptr;

struct pbc_param_s {
  pbc_param_interface_ptr api;
  void *data;
};
typedef struct pbc_param_s *pbc_param_ptr;
typedef struct pbc_param_s pbc_param_t[1];

static inline void pairing_init_pbc_param(struct pairing_s *pairing, pbc_param_ptr par) {
  par->api->init_pairing(pairing, par->data);
}

static inline void pbc_param_out_str(FILE *stream, pbc_param_ptr par) {
  par->api->out_str(stream, par->data);
}

static inline void pbc_param_clear(pbc_param_ptr par) {
  par->api->clear(par->data);
}

#endif //__PBC_PARAM_H__
