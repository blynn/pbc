/* Implementation of Boneh-Waters broadcast encryption scheme
   Code by:  Matt Steiner   MattS@cs.stanford.edu
   bce.h
*/

#include "pbc.h"
#include <gmp.h>
#include <string.h>


/* **********************************************************
   DEBUG having the debug flag turned on spews out lots of
   debugging output.
*********************************************************  */
#define DEBUG 1


/* **********************************************************
   PRIVATE KEY STRUCT
   This struct stores a single users' group elements and 
   their private key.  It also contains the recipients 
   currently in their product (a bit-vector representation),
   their decryption product (excluding their group element),
   and their index number.
********************************************************** */
typedef struct single_priv_key_s {
  element_t g_i_gamma;
  element_t g_i;
  element_t h_i;
  element_t decr_prod;
  int index;
}* priv_key_t;


/* **********************************************************
   GLOBAL BROADCAST PARAMS--
   Stores the:
   curve info--PUBLIC
   group elements--PUBLIC
   num-users-PUBLIC
*********************************************************  */
typedef struct global_broadcast_params_s {
  pairing_t pairing;
  char *pairFileName;
  element_t g;
  element_t h;
  element_t *gs;
  element_t *hs;
  int num_users;
}* global_broadcast_params_t;

/* **********************************************************
   BROADCAST SYSTEM stores:
   encryption product - can be public
   public key - public
   priv key - private
*********************************************************  */
typedef struct broadcast_system_s {
  element_t encr_prod;
  element_t pub_key;
  mpz_t priv_key;
}* broadcast_system_t;

/* **********************************************************
   CIPHERTEXT STRUCT
   Contains two group elements HDR C0 and HDR C1
*********************************************************  */
typedef struct ciphertext_s {
  element_t C0;
  element_t C1;
}* ct_t;

/* **********************************************************
   These functions free the memory associated with various 
   structures.  Note that the pointer you pass in will not
   be freed--you must free it manually to prevent freeing
   stack memory.
********************************************************** */

void FreeCT(ct_t myct);
void FreeBCS(broadcast_system_t bcs);
void FreeGBP(global_broadcast_params_t gbp);
void FreePK(priv_key_t key);


/* **********************************************************
   Sets up a global broadcast system by generating all of
   the gs, the hs, and their inverses.  Chooses random alpha
   for the exponent.  num_users must be a multiple of 8.
*********************************************************  */
void Setup_global_broadcast_params(global_broadcast_params_t *gbp, 
				   int num_users, char *pairFileName);

/* **********************************************************
   Stores the global broadcast system parameters to a file.
   WARNING: FILE WILL BE LARGE for large numbers of users
*********************************************************  */
void StoreParams(char *systemFileName, 
		 global_broadcast_params_t gbp,
		 broadcast_system_t sys);


/* **********************************************************
   Loads the global broadcast system paramters from a file.
*********************************************************  */
void LoadParams(char *systemFileName,
		global_broadcast_params_t *gbp,
		broadcast_system_t *sys);


/* **********************************************************
   Stores a single private key to a file.  The pairing file
   should be distributed with the private key file.
**********************************************************  */
void StorePrivKey(char *keyFileName, priv_key_t mykey);


/* **********************************************************
   Loads a single private key into a private key structure.
   Should be done after loading the pairing file.
**********************************************************  */
void LoadPrivKey(char *keyFileName, priv_key_t *mykey,
		 global_broadcast_params_t gbp);


/* **********************************************************
   Sets up one instance of a broadcast system, generating
   system specific global public and private keys.
*********************************************************  */
void Gen_broadcast_system(global_broadcast_params_t gbp, 
			  broadcast_system_t *sys);


/* **********************************************************
   This function gets the private key for a user at index i.
*********************************************************  */
void Get_priv_key(global_broadcast_params_t gbp, 
		  broadcast_system_t sys,
		  int i, priv_key_t mykey);


/* **********************************************************
   This function generates the encryption product from
   a bit vector representing the users who should be able
   to decrypt the message.
*********************************************************  */
void Gen_encr_prod_from_bitvec(global_broadcast_params_t gbp, 
			       broadcast_system_t sys, char *recip);


/* **********************************************************
   This function generates the encryption product from
   an array of indicies corresponding to the users who 
   should be able to decrypt the message.  You must
   pass in the correct array size.
*********************************************************  */
void Gen_encr_prod_from_indicies(global_broadcast_params_t gbp, 
				 broadcast_system_t sys,
				 int *in_recip, int num_recip);


/* **********************************************************
   This function changes the encryption product of the system
   first by removing the N_rems elements in the rems array
   and then by adding the N_adds elements in the adds array.
*********************************************************  */
void Change_encr_prod_indicies(global_broadcast_params_t gbp, 
			       broadcast_system_t sys,
			       int *adds, int N_adds, 
			       int *rems, int N_rems);


/* **********************************************************
   This function generates the decryption product from
   a bit vector representing the users who should be able
   to decrypt the message.  The product gets stored into 
   mykey.  Group element corresponding to receiver is not
   included.
*********************************************************  */
void Gen_decr_prod_from_bitvec(global_broadcast_params_t gbp, 
			       int receiver, 
			       char *recip, priv_key_t mykey);



/* **********************************************************
   This function generates the decryption product from
   an array of indicies corresponding to the users who 
   should be able to decrypt the message.  You must
   pass in the correct array size.  The product gets stored 
   into mykey.  Group element corresponding to receiver is 
   not included.
*********************************************************  */
void Gen_decr_prod_from_indicies(global_broadcast_params_t gbp, 
				 int receiver, int *in_recip, 
				 int num_recip, priv_key_t mykey);


/* **********************************************************
   This function changes the decryption product of the user
   first by removing the N_rems elements in the rems array
   and then by adding the N_adds elements in the adds array.
   You must pass in the correct array size.  The product gets 
   stored into mykey.  Group element corresponding to receiver 
   is not included.
*********************************************************  */
void Change_decr_prod_indicies(global_broadcast_params_t gbp, 
			       int receiver, 
			       int *adds, int N_adds, 
			       int *rems, int N_rems, priv_key_t mykey);


/* **********************************************************
   An extremely useful function for validating your results
   This function takes a bit-string and prints out the first
   length bytes of it.  Accordingly, you must give the
   function length = Num_Users/8
*********************************************************  */
void PrintBitString(char *bs, int length);


/* **********************************************************
   This function generates a broadcast key and a cipher-text
   header, once the encryption product has been calculated.
*********************************************************  */
void BroadcastKEM_using_product(global_broadcast_params_t gbp, 
				broadcast_system_t sys,
				ct_t myct, element_t key);


/* **********************************************************
   This function generates a broadcast key and a cipher-text
   header, by calling the Gen-prod-from-bitvec and then
   BroadcastKEM-using-product routines.  Just a wrapper.
*********************************************************  */
void BroadcastKEM_using_bitvec(global_broadcast_params_t gbp,
			       broadcast_system_t sys, 
			       char *recip, ct_t myct, element_t key);


/* **********************************************************
   This function generates a broadcast key and a cipher-text
   header, by calling the Gen-prod-from-indicies and then
   BroadcastKEM-using-product routines.  Just a wrapper.
*********************************************************  */
void BroadcastKEM_using_indicies(global_broadcast_params_t gbp,
				 broadcast_system_t sys, ct_t myct,
				 int *in_recip, int num_recip, 
				 element_t key);


/* **********************************************************
   This function retrieves a broadcast key from a cipher-text
   header, once the decryption product has been calculated.
*********************************************************  */
void DecryptKEM_using_product(global_broadcast_params_t gbp, 
			      priv_key_t mykey, element_t key,
			      ct_t myct);



/* **********************************************************
   This function retrieves a broadcast key from a cipher-text
   header, once the decryption product has been calculated, 
   by calling the Gen-decr-prod-from-bitvec and then
   DecryptKEM-using-product routines.  Just a wrapper.
*********************************************************  */
void Decrypt_BC_KEM_using_bitvec(global_broadcast_params_t gbp,
				 priv_key_t mykey, element_t key,
				 char *recip, ct_t myct);


/* **********************************************************
   This function retrieves a broadcast key from a cipher-text
   header, once the decryption product has been calculated, 
   by calling the Gen-decr-prod-from-bitvec and then
   DecryptKEM-using-product routines.  Just a wrapper.
*********************************************************  */
void Decrypt_BC_KEM_using_indicies(global_broadcast_params_t gbp, 
				   priv_key_t mykey, element_t key,
				   ct_t myct, int *in_recip, 
				   int num_recip);


