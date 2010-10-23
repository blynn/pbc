// Requires:
// * stdio.h
// * gmp.h
// * utils.h
// * field.h
// * param.h
#ifndef __PBC_PAIRING_H__
#define __PBC_PAIRING_H__

struct pairing_pp_s {
  struct pairing_s *pairing;
  void *data;
};
typedef struct pairing_pp_s pairing_pp_t[1];
typedef struct pairing_pp_s *pairing_pp_ptr;

struct pairing_s {
  mpz_t r; // order of G1, G2, GT
  field_t Zr; // the field Z_r
  field_ptr G1, G2;
  field_t GT; // group of rth roots of unity

  mpz_t phikonr;
  // Phi_k(q)/r where Phi_k is the kth cyclotomic polynomial,
  // q as in F_q, is the base field

  void (*phi)(element_ptr out, element_ptr in, struct pairing_s *pairing); //isomorphism G2 --> G1
  void (*map)(element_ptr out, element_ptr in1, element_ptr in2,
      struct pairing_s *p);
  void (*prod_pairings)(element_ptr out, element_t in1[], element_t in2[], int n_prod,
            struct pairing_s *p);  //calculate a product of pairings at one time.
  // is_almost coddh returns true given (g, g^x, h, h^x) or (g, g^x, h, h^-x)
  // order is important: a, b are from G1, c, d are from G2
  int (*is_almost_coddh)(element_ptr a, element_ptr b,
      element_ptr c, element_ptr d,
      struct pairing_s *p);
  void (*clear_func)(struct pairing_s *);
  void (*pp_init)(pairing_pp_t p, element_t in1, struct pairing_s *);
  void (*pp_clear)(pairing_pp_t p);
  void (*pp_apply)(element_t out, element_t in2, pairing_pp_t p);
  void (*finalpow)(element_t e);
  void (*option_set)(struct pairing_s *, char *key, char *value);
  void *data;
};

typedef struct pairing_s pairing_t[1];
typedef struct pairing_s *pairing_ptr;

// TODO: The 'pairing' argument is redundant.
/*@manual pairing_apply
Get ready to perform a pairing whose first input is 'in1',
and store the results of time-saving precomputation in 'p'.
*/
static inline void pairing_pp_init(pairing_pp_t p, element_t in1, pairing_t pairing) {
  if (element_is0(in1)) {
    p->pairing = NULL;
    return;
  }
  p->pairing = pairing;
  pairing->pp_init(p, in1, pairing);
}

/*@manual pairing_apply
Clear 'p'. This should be called after 'p' is no longer needed.
*/
static inline void pairing_pp_clear(pairing_pp_t p) {
  if (!p->pairing) {
    // happens when p was initialized with identity
    return;
  }
  p->pairing->pp_clear(p);
}

/*@manual pairing_apply
Compute a pairing using 'in2' and the preprocessed information stored in 'p'
and store the output in 'out'. The inputs to the pairing are the element
previously used to initialize 'p' and the element 'in2'.
*/
static inline void pairing_pp_apply(element_t out, element_t in2, pairing_pp_t p) {
  if (!p->pairing) {
    element_set0(out);
    return;
  }
  if (element_is0(in2)) {
    element_set0(out);
    return;
  }
  p->pairing->pp_apply((element_ptr) out->data, in2, p);
}

/*@manual pairing_init
Initialize pairing from parameters in a ASCIIZ string 'str'
Returns 0 on success, 1 on failure.
*/
int pairing_init_set_str(pairing_t pairing, const char *s);

/*@manual pairing_init
Same, but read at most 'len' bytes.
If 'len' is 0, it behaves as the previous function.
Returns 0 on success, 1 on failure.
*/
int pairing_init_set_buf(pairing_t pairing, const char *s, size_t len);

/*@manual pairing_init
Initialize a pairing with pairing parameters 'p'.
*/
void pairing_init_pbc_param(struct pairing_s *pairing, pbc_param_ptr p);

/*@manual pairing_init
Free the space occupied by 'pairing'. Call
whenever a +pairing_t+ variable is no longer needed.
Only call this after all elements associated with 'pairing'
have been cleared, as they need information stored in the 'pairing'
structure.
*/
void pairing_clear(pairing_t pairing);

static inline void pairing_apply(element_t out, element_t in1, element_t in2,
    pairing_t pairing) {
  PBC_ASSERT(pairing->GT == out->field, "pairing output mismatch");
  PBC_ASSERT(pairing->G1 == in1->field, "pairing 1st input mismatch");
  PBC_ASSERT(pairing->G2 == in2->field, "pairing 2nd input mismatch");
  if (element_is0(in1)) {
    element_set0(out);
    return;
  }
  if (element_is0(in2)) {
    element_set0(out);
    return;
  }
  // TODO: 'out' is an element of a multiplicative subgroup, but the
  // pairing routine expects it to be an element of the full group, hence
  // the 'out->data'. I should make this clearer.
  pairing->map((element_ptr) out->data, in1, in2, pairing);
}

/*@manual pairing_apply
Computes a pairing: 'out' = 'e'('in1', 'in2'),
where 'in1', 'in2', 'out' must be in the groups G1, G2, GT.
*/
static inline void element_pairing(element_t out, element_t in1, element_t in2) {
  pairing_ptr pairing = out->field->pairing;
  PBC_ASSERT(pairing != NULL, "pairing output mismatch");
  pairing_apply(out, in1, in2, pairing);
}

/*@manual pairing_apply
Computes the product of pairings, that is
'out' = 'e'('in1'[0], 'in2'[0]) ... 'e'('in1'[n-1], 'in2'[n-1]).
The arrays 'in1', 'in2' must have at least 'n' elements belonging to
the groups G1, G2 respectively, and 'out' must belong to the group GT.
*/
static inline void element_prod_pairing(
    element_t out, element_t in1[], element_t in2[], int n) {
  pairing_ptr pairing = out->field->pairing;
  int i;
  PBC_ASSERT(pairing->GT == out->field, "pairing output mismatch");
  for(i = 0; i < n; i++) {
    PBC_ASSERT(pairing->G1 == in1[i]->field, "pairing 1st input mismatch");
    PBC_ASSERT(pairing->G2 == in2[i]->field, "pairing 2nd input mismatch");
    if (element_is0(in1[i])) {
      element_set0(out);
      return;
    }
    if (element_is0(in2[i])) {
      element_set0(out);
      return;
    }
  }
  pairing->prod_pairings((element_ptr) out->data, in1, in2, n, pairing);
}

/*@manual pairing_op
Returns true if G1 and G2 are the same group.
*/
static inline int pairing_is_symmetric(pairing_t pairing) {
  return pairing->G1 == pairing->G2;
}

/*@manual pairing_op
Returns the length in bytes needed to represent an element of G1.
*/
static inline int pairing_length_in_bytes_G1(pairing_t pairing) {
  return pairing->G1->fixed_length_in_bytes;
}

/*@manual pairing_op
Returns the length in bytes needed to represent the x-coordinate of
an element of G1.
*/
static inline int pairing_length_in_bytes_x_only_G1(pairing_t pairing) {
  return pairing->G1->fixed_length_in_bytes / 2;
}

/*@manual pairing_op
Returns the length in bytes needed to represent a compressed form of
an element of G1. There is some overhead in decompressing.
*/
static inline int pairing_length_in_bytes_compressed_G1(pairing_t pairing) {
  return pairing->G1->fixed_length_in_bytes / 2 + 1;
}

/*@manual pairing_op
Returns the length in bytes needed to represent an element of G2.
*/
static inline int pairing_length_in_bytes_G2(pairing_t pairing) {
  return pairing->G2->fixed_length_in_bytes;
}

/*@manual pairing_op
Returns the length in bytes needed to represent a compressed form of
an element of G2. There is some overhead in decompressing.
*/
static inline int pairing_length_in_bytes_compressed_G2(pairing_t pairing) {
  return pairing->G2->fixed_length_in_bytes / 2 + 1;
}

/*@manual pairing_op
Returns the length in bytes needed to represent the x-coordinate of
an element of G2.
*/
static inline int pairing_length_in_bytes_x_only_G2(pairing_t pairing) {
  return pairing->G2->fixed_length_in_bytes / 2;
}

/*@manual pairing_op
Returns the length in bytes needed to represent an element of GT.
*/
static inline int pairing_length_in_bytes_GT(pairing_t pairing) {
  return pairing->GT->fixed_length_in_bytes;
}

/*@manual pairing_op
Returns the length in bytes needed to represent an element of Zr.
*/
static inline int pairing_length_in_bytes_Zr(pairing_t pairing) {
  return pairing->Zr->fixed_length_in_bytes;
}

static inline int is_almost_coddh(element_t a, element_t b,
    element_t c, element_t d, pairing_t pairing) {
  return pairing->is_almost_coddh(a, b, c, d, pairing);
}

/*@manual einit.1
*/
static inline void element_init_G1(element_t e, pairing_t pairing) {
  element_init(e, pairing->G1);
}

/*@manual einit.1
*/
static inline void element_init_G2(element_t e, pairing_t pairing) {
  element_init(e, pairing->G2);
}

/*@manual einit.1
Initialize 'e' to be an element of the group G1, G2 or GT of 'pairing'.
*/
static inline void element_init_GT(element_t e, pairing_t pairing) {
  element_init(e, pairing->GT);
}

/*@manual einit.1
Initialize 'e' to be an element of the ring Z_r of 'pairing'.
r is the order of the groups G1, G2 and GT that are involved in the pairing.
*/
static inline void element_init_Zr(element_t e, pairing_t pairing) {
  element_init(e, pairing->Zr);
}

static inline void pairing_option_set(pairing_t pairing, char *key, char *value) {
  pairing->option_set(pairing, key, value);
}

// Initialize GT = group of rth roots of unity in f.
// Requires pairing->r has been set.
void pairing_GT_init(pairing_ptr pairing, field_t f);

#endif //__PBC_PAIRING_H__
