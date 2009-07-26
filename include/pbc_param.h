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

/*@manual param
Initialize a pairing with pairing parameters 'p'.
*/
static inline void pairing_init_pbc_param(struct pairing_s *pairing, pbc_param_ptr p) {
  p->api->init_pairing(pairing, p->data);
}

/*@manual param
Write pairing parameters to ''stream'' in a text format.
*/
static inline void pbc_param_out_str(FILE *stream, pbc_param_ptr p) {
  p->api->out_str(stream, p->data);
}

/*@manual param
Clear 'p'. Call after 'p' is no longer needed.
*/
static inline void pbc_param_clear(pbc_param_ptr p) {
  p->api->clear(p->data);
}

#endif //__PBC_PARAM_H__
