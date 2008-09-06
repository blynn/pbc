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
For a given discriminant D, searches for MNT curves of embedding degree 6
suitable for cryptography (type D pairings) where the group order
is at most 'bitlimit' bits. For each suitable set of CM parameters found,
call supplied callback with +cm_info_t+ and given void pointer.
Returns the number of CM parameters found.
*/
int find_mnt6_curve(void (*callback)(cm_info_ptr, void *), void *data,
    unsigned int D, unsigned int bitlimit);

/*@manual cminfo
For a given discriminant D, searches for a Freeman curve of embedding
degree 10
suitable for cryptography (type D pairings) where the group order
is at most 'bitlimit' bits. For each suitable set of CM parameters found,
call supplied callback with +cm_info_t+ and given void pointer.
Returns the number of CM parameters found.
*/
int find_freeman_curve(void (*callback)(cm_info_ptr, void *), void *data,
    unsigned int D, unsigned int bitlimit);

#endif //__PBC_MNT_H__
