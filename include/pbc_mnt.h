//requires
// * gmp.h
#ifndef __PBC_MNT_H__
#define __PBC_MNT_H__

struct cm_info_s {
  mpz_t q; //curve defined over F_q
  mpz_t n; //has order n (= q - t + 1) in F_q (and r^2 in F_q^k)
  mpz_t h; //h * r = n, r is prime
  mpz_t r;
  int D; //discrminant needed to find j-invariant
  int k; //embedding degree
};

typedef struct cm_info_s *cm_info_ptr;
typedef struct cm_info_s cm_info_t[1];

/*@manual cminfo
Initializes 'cm'.
*/
void cm_info_init(cm_info_t cm);
/*@manual cminfo
Clears 'cm'.
*/
void cm_info_clear(cm_info_t cm);

/*@manual cminfo
For a given discriminant D, searches for type D pairings suitable for
cryptography (MNT curves of embedding degree 6).
The group order is at most 'bitlimit' bits. For each set of CM parameters
found, call 'callback' with +cm_info_t+ and given +void *+. If the callback
returns nonzero, stops search and returns that value.
Otherwise returns 0.
*/
int cm_search_d(int (*callback)(cm_info_ptr, void *), void *data,
  unsigned int D, unsigned int bitlimit);

/*@manual cminfo
For a given discriminant D, searches for type G pairings suitable for
cryptography (Freeman curve).
The group order is at most 'bitlimit' bits. For each set of CM parameters
found, call 'callback' with +cm_info_t+ and given +void *+. If the callback
returns nonzero, stops search and returns that value.
Otherwise returns 0.
*/
int cm_search_g(int (*callback)(cm_info_ptr, void *), void *data,
  unsigned int D, unsigned int bitlimit);

#endif //__PBC_MNT_H__
