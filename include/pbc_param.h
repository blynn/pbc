// Requires:
// * gmp.h
#ifndef __PBC_PARAM_H__
#define __PBC_PARAM_H__

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
