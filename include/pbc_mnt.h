//requires
// * gmp.h
// * darray.h
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
Initializes ''cm''.
*/
void cm_info_init(cm_info_t cm);
/*@manual cminfo
Clears ''cm''.
*/
void cm_info_clear(cm_info_t cm);

/*@manual cminfo
For a given discriminant D, searches for MNT curves of embedding degree 6
suitable for cryptography (type D pairings) where the group order
is at most ''bitlimit'' bits. For each suitable set of CM parameters found,
a <type>cm_info_t</type> is created and appended to the dynamic array ''L''.
Returns the number of CM parameters found.
</para>
<para>
When no longer needed, <command>cm_info_clear</command> should be
called on each appended element in ''L''. (And of course, when ''L'' is
no longer needed <command>darray_clear</command> should be called on ''L''.)
*/
int find_mnt6_curve(darray_t L, unsigned int D, unsigned int bitlimit);

/*@manual cminfo
For a given discriminant D, searches for a Freeman curve of embedding
degree 10
suitable for cryptography (type D pairings) where the group order
is at most ''bitlimit'' bits. For each suitable set of CM parameters found,
a <type>cm_info_t</type> is created and appended to the dynamic array ''L''.
Returns the number of CM parameters found.
</para>
<para>
When no longer needed, <command>cm_info_clear</command> should be
called on each appended element in ''L''. (And of course, when ''L'' is
no longer needed <command>darray_clear</command> should be called on ''L''.)
*/
int find_freeman_curve(darray_t L, unsigned int D, unsigned int bitlimit);

#endif //__PBC_MNT_H__
