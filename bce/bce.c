/* Implementation of Boneh-Waters broadcast encryption scheme
 * Code by:  Matt Steiner   MattS@cs.stanford.edu
 * bce.c
 */

#include "pbc.h"
#include <gmp.h>
#include <string.h>
#include <stdio.h>
#include "bce.h"
#include <stdlib.h>


void FreeCT(ct_t myct)
{
  if(!myct) {
    printf("error: null pointer passed to freeCT\n");
    return;
  }
  element_clear(myct->C0);
  element_clear(myct->C1);
  return;
}

void FreeBCS(broadcast_system_t bcs)
{
  if(!bcs) {
    printf("error: null pointer passed to freeBCS\n");
    return;
  }
  element_clear(bcs->encr_prod);
  element_clear(bcs->pub_key);  
  mpz_clear(bcs->priv_key);
  return;
}

void FreeGBP(global_broadcast_params_t gbp)
{
  if(!gbp) {
    printf("error: null pointer passed to freeGBP\n");
    return;
  }
  //do something about the pairing
  free(gbp->pairFileName);
  element_clear(gbp->g);
  element_clear(gbp->h);
  int i;
  for(i = 0; i < gbp->num_users; i++) {
    if(i == gbp->num_users) continue;
    element_clear(gbp->gs[i]);
    element_clear(gbp->hs[i]);
  }
  free(gbp->gs);
  free(gbp->hs);
  
  return;
}

void FreePK(priv_key_t key)
{
  if(!key) {
    printf("error: null pointer passed to freePK\n");
    return;
  }
  element_clear(key->g_i_gamma);
  element_clear(key->g_i);
  element_clear(key->h_i);
  element_clear(key->decr_prod);
  return;
}


static inline void out(element_t elem, FILE *myfile) 
{
  int sz = element_length_in_bytes(elem);
  fwrite(&sz, 4, 1, myfile);
  unsigned char* data = malloc(sz);
  if(!data) printf("DATA IS NULL\n");
  element_to_bytes(data, elem);
  fwrite(data, sz, 1, myfile);
  free(data);
}

static inline void in(element_t elem, FILE *myfile) {
  int sz;
  fread(&sz, 4, 1, myfile);
  unsigned char* data = malloc(sz);
  fread(data, sz, 1, myfile);
  element_from_bytes(elem, data);
  free(data);
}

void StorePrivKey(char *keyFileName, priv_key_t mykey)
{
  if(!mykey) {
    printf("ACK!  You gave me no key!  I die.\n");
    return;
  }
  if(!keyFileName) {
    printf("ACK!  You gave me no key filename!  I die.\n");
    return;
  }
  FILE *keyf = fopen(keyFileName, "w");
  if(!keyf) {
    printf("ACK! couldn't write to file system.  I die\n");
    return;
  }
  //store key index
  fwrite(&(mykey->index),4,1, keyf);
  //if(DEBUG) printf("done storing key index\n");
  
  //store g_i_gamma
  out(mykey->g_i_gamma, keyf);
  //if(DEBUG) printf("done storing g_i_gamma\n");

  //store g_i
  out(mykey->g_i, keyf);
  //if(DEBUG) printf("done storing g_i\n");

  //store h_i
  out(mykey->h_i, keyf);
  //if(DEBUG) printf("done storing h_i\n");

  //store decr_prod
  out(mykey->decr_prod, keyf);
  //if(DEBUG) printf("done storing decr_prod\n");

  fclose(keyf);
  return;
}


void LoadPrivKey(char *keyFileName, priv_key_t *mykey, 
		 global_broadcast_params_t gbp)
{
  if(!gbp) {
    printf("ACK!  You gave me no broadcast params!  I die.\n");
    return;
  }  
  if(!mykey) {
    printf("ACK!  You gave me no key!  I die.\n");
    return;
  }
  if(!keyFileName) {
    printf("ACK!  You gave me no key filename!  I die.\n");
    return;
  }
  FILE *keyf = fopen(keyFileName, "r");
  if(!keyf) {
    printf("ACK! couldn't write to file system.  I die\n");
    return;
  }
  
  priv_key_t key = (priv_key_t) malloc(sizeof(struct single_priv_key_s));
  
  //restore key index
  fread(&(key->index),4,1, keyf);
  //if(DEBUG) printf("done restoring key index\n");
  
  //restore g_i_gamma
  element_init(key->g_i_gamma, gbp->pairing->G1);
  in(key->g_i_gamma, keyf);
  //if(DEBUG) printf("done restoring g_i_gamma\n");

  //restore g_i
  element_init(key->g_i, gbp->pairing->G1);
  in(key->g_i, keyf);
  //if(DEBUG) printf("done restoring g_i\n");

  //restore h_i
  element_init(key->h_i, gbp->pairing->G2);
  in(key->h_i, keyf);
  //if(DEBUG) printf("done restoring h_i\n");

  //restore decr_prod
  element_init(key->decr_prod, gbp->pairing->G1);
  in(key->decr_prod, keyf);
  //if(DEBUG) printf("done restoring decr_prod\n");

  fclose(keyf);
  *mykey = key;
  return;

}



void StoreParams(char *systemFileName, 
		 global_broadcast_params_t gbp,
		 broadcast_system_t sys)
{

  if(!gbp) {
    printf("ACK!  You gave me no broadcast params!  I die.\n");
    return;
  }
  if(!sys) {
    printf("ACK!  You gave me no broadcast system!  I die.\n");
    return;
  }
  if(!systemFileName) {
    printf("ACK!  You gave me no system filename!  I die.\n");
    return;
  }

  FILE *sysp = fopen(systemFileName, "w");
  if(!sysp) {
    printf("ACK! couldn't write to file system.  I die\n");
    return;
  }

  int leng = strlen(gbp->pairFileName) + 1;
  fwrite(&leng, 4, 1, sysp);
  fwrite(gbp->pairFileName, 1, leng, sysp);

  //store num_users
  fwrite(&(gbp->num_users),4,1, sysp);
  //if(DEBUG) printf("done storing n_users\n");

  //store encr_prod
  out(sys->encr_prod, sysp);
  //if(DEBUG) printf("done storing encr_prod\n");

  //store pub_key
  out(sys->pub_key, sysp);
  //if(DEBUG) printf("done storing pub_key\n");

  //store g 
  out(gbp->g, sysp);
  //if(DEBUG) printf("done storing g\n");

  //store gs
  int i;
  for(i = 0; i < 2*gbp->num_users; i++) {
    if(i == gbp->num_users) continue;
    out(gbp->gs[i], sysp);
    //if(DEBUG) printf("done storing g %d\n",i);    
  }
  //if(DEBUG) printf("done storing gs\n");

  //store h
  out(gbp->h, sysp);
  //if(DEBUG) printf("done storing h\n");

  //store hs
  for(i = 0; i < 2*gbp->num_users; i++) {
    if(i == gbp->num_users) continue;
    out(gbp->hs[i], sysp);
    //if(DEBUG) printf("done storing h %d\n",i);
  }
  //if(DEBUG) printf("done storing hs\n");

  fclose(sysp);

  return;

}


void LoadParams(char *systemFileName,
		global_broadcast_params_t *gbp,
		broadcast_system_t *sys)
{

  if(!gbp) {
    printf("ACK!  You gave me no broadcast params!  I die.\n");
    return;
  }
  if(!gbp) {
    printf("ACK!  You gave me no broadcast system!  I die.\n");
    return;
  }
  if(!systemFileName) {
    printf("ACK!  You gave me no system filename!  I die.\n");
    return;
  }

  global_broadcast_params_t p;
  broadcast_system_t s;

  p = malloc(sizeof(struct global_broadcast_params_s));
  s = malloc(sizeof(struct broadcast_system_s));

  FILE *sysp = fopen(systemFileName, "r");
  if(!sysp) {
    printf("ACK! couldn't open %s  I die\n", systemFileName);
    return;
  }

  int leng;
  fread(&leng, 4, 1, sysp);
  p->pairFileName = (char *) malloc(leng);
  fread(p->pairFileName, 1, leng, sysp);
  FILE *params = fopen(p->pairFileName, "r");
  if(!params) {
    printf("ACK! couldn't open %s  I die\n", p->pairFileName);
    return;    
  }
  pairing_init_inp_str(p->pairing, params);
  fclose(params);

  //restore num_users
  fread(&(p->num_users),4,1, sysp); 

  //restore encr_prod
  element_init(s->encr_prod, p->pairing->G1);
  in(s->encr_prod, sysp);
  //element_out_str(stdout, 0, s->encr_prod);

  //restore pub_key
  element_init(s->pub_key, p->pairing->G1);
  in(s->pub_key, sysp);
  //element_out_str(stdout, 0, s->pub_key);

  //restore g 
  element_init(p->g, p->pairing->G1);
  in(p->g, sysp);


  p->gs = malloc(2 * p->num_users * sizeof(element_t));
  p->hs = malloc(2 * p->num_users * sizeof(element_t));

  //restore gs
  int i;
  for(i = 0; i < 2*p->num_users; i++) {
    if(i == p->num_users) continue;
    element_init(p->gs[i], p->pairing->G1);
    in(p->gs[i], sysp);
  }
  
  //restore h
  element_init(p->h, p->pairing->G2);
  in(p->h, sysp);
  
  //restore hs
  for(i = 0; i < 2*p->num_users; i++) {
    if(i == p->num_users) continue;
    element_init(p->hs[i], p->pairing->G2);
    in(p->hs[i], sysp);
  }
    
  fclose(sysp);

  //now insert a dummy private key
  mpz_init(s->priv_key);
  
  *gbp = p;
  *sys = s;

  return;


}

void DecryptKEM_using_product(global_broadcast_params_t gbp, 
			      priv_key_t mykey, element_t key,
			      ct_t myct)
{
  if(!gbp) {
    printf("ACK!  You gave me no broadcast params!  I die.\n");
    return;
  }
  if(!mykey) {
    printf("ACK!  You gave me no private key info  I die.\n");
    return;
  }
  if(!myct) {
    printf("ACK!  No struct cipher text to decode  I die.\n");
    return;
  }
  if(!key) {
    printf("ACK!  No place to put my key!  I die.\n");
    return;
  }
  if(!mykey->decr_prod) {
     printf("ACK!  Calculate decryption prodcut before ");
     printf("calling this function! I die.\n");
     return;
  }
  
  element_t temp;
  element_t temp2;
  element_t di_de;
  element_t temp3;
  
  element_init(temp, gbp->pairing->GT);
  element_init(temp2, gbp->pairing->GT);
  element_init(di_de, gbp->pairing->G1);
  element_init(temp3, gbp->pairing->GT);
 
  //Generate the numerator
  bilinear_map(temp, myct->C1, mykey->h_i, gbp->pairing);
  //G1 element in denom
  element_mul(di_de, mykey->g_i_gamma, mykey->decr_prod);
  //Generate the denominator
  bilinear_map(temp2, di_de, myct->C0, gbp->pairing);
  //Invert the denominator
  element_invert(temp3, temp2);

  
  element_init(key, gbp->pairing->GT);
  //multiply the numerator by the inverted denominator
  element_mul(key, temp, temp3);
  
}

void Decrypt_BC_KEM_using_bitvect(global_broadcast_params_t gbp, 
				  priv_key_t mykey, element_t key,
				  ct_t myct, char *recip)
{
  if(!mykey) {
    printf("\nyou didn't give me a valid key.  I die\n");
    return;
  }
  Gen_decr_prod_from_bitvec(gbp, mykey->index, recip, mykey);
  DecryptKEM_using_product(gbp, mykey, key, myct);
}


void Decrypt_BC_KEM_using_indicies(global_broadcast_params_t gbp, 
				   priv_key_t mykey, element_t key,
				   ct_t myct, int *in_recip, int num_recip)
{
  if(!mykey) {
    printf("\nyou didn't give me a valid key.  I die\n");
    return;
  } 
  Gen_decr_prod_from_indicies(gbp, mykey->index,in_recip, num_recip, mykey);
  DecryptKEM_using_product(gbp, mykey, key, myct);  
}

			      


void BroadcastKEM_using_bitvec(global_broadcast_params_t gbp,
			       broadcast_system_t sys, 
			       char *recip, ct_t myct, element_t key)
{

  Gen_encr_prod_from_bitvec(gbp, sys, recip);
  if(DEBUG && 0) {
    printf("bitvec product = ");
    element_out_str(stdout, 0, sys->encr_prod);
    printf("\n");
  }
  BroadcastKEM_using_product(gbp, sys, myct, key);
}

void BroadcastKEM_using_indicies(global_broadcast_params_t gbp,
				 broadcast_system_t sys, ct_t myct,
				 int *in_recip, int num_recip, element_t key)
{
  Gen_encr_prod_from_indicies(gbp, sys, in_recip, num_recip);
  if(DEBUG && 0) {
    printf("index product = ");
    element_out_str(stdout, 0, sys->encr_prod);
    printf("\n");
  }
  BroadcastKEM_using_product(gbp, sys, myct, key);
}




void BroadcastKEM_using_product(global_broadcast_params_t gbp, 
				broadcast_system_t sys,
				ct_t myct, element_t key)
{

  if(!gbp) {
    printf("ACK!  You gave me no broadcast params!  I die.\n");
    return;
  }
  if(!sys) {
    printf("ACK!  You gave me no broadcast system!  I die.\n");
    return;
  }
  if(!myct) {
    printf("ACK!  No struct to store return vals!  I die.\n");
    return;
  }

  mpz_t t;
  mpz_init(t);
  pbc_mpz_random(t, gbp->pairing->r);
  if(mpz_sgn(t) < 0) {
    mpz_neg(t,t);
  }
  
  element_init(key, gbp->pairing->GT);
  element_init(myct->C0, gbp->pairing->G2);
  element_init(myct->C1, gbp->pairing->G1);
  
  //COMPUTE K
  bilinear_map(key, gbp->gs[gbp->num_users-1], gbp->hs[0], gbp->pairing);
  element_pow(key, key, t);

  //COMPUTE C0
  element_pow(myct->C0, gbp->h, t);

  //COMPUTE C1
  if(DEBUG && 0) {
    printf("\npub_key = ");
    element_out_str(stdout, 0, sys->pub_key);
    printf("\nencr_prod = ");
    element_out_str(stdout, 0, sys->encr_prod);
  }
  element_mul(myct->C1, sys->pub_key, sys->encr_prod);
  if(DEBUG && 0) {
    printf("\npub_key = ");
    element_out_str(stdout, 0, sys->pub_key);
    printf("\nencr_prod = ");
    element_out_str(stdout, 0, sys->encr_prod);
    printf("\nhdr_c1 = ");
    element_out_str(stdout, 0, myct->C1);    
    printf("\n");
  }
  element_pow(myct->C1, myct->C1, t);


}

void Change_decr_prod_indicies(global_broadcast_params_t gbp, int receiver, 
			       int *adds, int N_adds, int *rems, int N_rems, 
			       priv_key_t mykey)
{
  // REMOVES THE OLD ONES, THEN ADDS THE NEW ONES
  int i;
  element_t temp_inv;
  int incl_num;
  if(!gbp) {
    printf("ACK!  You gave me no broadcast params!  I die.\n");
    return;
  }
  if(!mykey) {
    printf("ACK!  You gave me no broadcast system!  I die.\n");
    return;
  }
  int n = gbp->num_users;
  element_init(temp_inv, gbp->pairing->G1);
  if(rems) {
    for(i = 0; i < N_rems; i++) {
      //removing elements from the set after 
      //checking if it's in the bit-vector
      incl_num = rems[i];
      if(incl_num < 1 || incl_num > gbp->num_users) {
	printf("element %d was outside the range of valid users\n",i);
	printf("only give me valid values.  i die.\n");
	return;
      }
      if(incl_num == receiver) {
	if(DEBUG) printf("incl_num == receiver, continuing\n");
	continue;
      }
      element_invert(temp_inv, gbp->gs[(n-incl_num)+receiver]);
      element_mul(mykey->decr_prod, mykey->decr_prod, temp_inv);
    } 
  }
  if(adds) {
    for(i = 0; i < N_adds; i++) {
      //adding elements to the set after checking if they're already there
      incl_num = adds[i];
      if(incl_num < 1 || incl_num > gbp->num_users) {
	printf("element %d was outside the range of valid users\n",i);
	printf("only give me valid values.  i die.\n");
	return;
      }
      if(incl_num == receiver) {
	if(DEBUG) printf("incl_num == receiver, continuing\n");
	continue;
      }
      element_mul(mykey->decr_prod, mykey->decr_prod, 
		  gbp->gs[(n-incl_num)+receiver]);
    }  
  }
}



void Gen_decr_prod_from_indicies(global_broadcast_params_t gbp, int receiver,
				 int *in_recip, int num_recip, 
				 priv_key_t mykey)
{
  if(!gbp) {
    printf("ACK!  You gave me no broadcast params!  I die.\n");
    return;
  }
  if(!mykey) {
    printf("ACK!  You gave me no private key to put decryption product in.\n");
    return;
  }
  if(!in_recip) {
    printf("ACK!  You gave me no recipient list!  I die.\n");
    return;
  }
  //RECIPIENT NUMBERS MUST BE IN NORMAL NOTATION e.g. 1 to N
  int i;
  element_init(mykey->decr_prod, gbp->pairing->G1);

  int n = gbp->num_users;

  int incl_num;
  int already_set = 0;

  for(i = 0; i < num_recip; i++) {
    incl_num = in_recip[i];
    if(incl_num < 1 || incl_num > gbp->num_users) {
      printf("element %d was outside the range of valid users\n",i);
      printf("only give me valid values.  i die.\n");
      return;
    }
    if(incl_num == receiver) {
      if(DEBUG && 0) printf("incl_num == receiver, continuing\n");
      continue;
    }
    if(DEBUG && 0) printf("\nputting element %d in product",incl_num);
    if(!already_set) {
      element_set(mykey->decr_prod, gbp->gs[(n-incl_num)+receiver]);
      already_set = 1;
    } else {
      element_mul(mykey->decr_prod, mykey->decr_prod, 
		  gbp->gs[(n-incl_num)+receiver]);
    }
  }
  
}


void Gen_decr_prod_from_bitvec(global_broadcast_params_t gbp, 
			       int receiver, char *recip, priv_key_t mykey)
{
  
  // recip must have length num_users/8;
  // BITVector[0    ]lsb corresponds to Recipient 1, and
  // BITVector[n/8-1]msb corresponds to Recipient n
  char working;
  int main_index = 1;
  int i,j;
  int already_set = 0;

  
  if(!gbp) {
    printf("ACK!  You gave me no broadcast params!  I die.\n");
    return;
  }
  if(!recip) {
    printf("ACK!  You gave me no recipient list!  I die.\n");
    return;
  }
  if(!mykey) {
    printf("ACK!  You gave me no private key to put decryption product in.\n");
  }
  element_init(mykey->decr_prod, gbp->pairing->G1);
  int n = gbp->num_users;
  int length = n / 8;
  
  for(i = 0; i < length; i++) {
    working = recip[i];
    for(j = 0; j < 8; j++) {
      if(main_index == receiver) {
	main_index++;
	working = working>>1;
	continue;
      }
      if(working & 1) {
	if(!already_set) {
	  element_set(mykey->decr_prod, gbp->gs[(n-main_index)+receiver]);
	  already_set = 1;
	} else {
	  element_mul(mykey->decr_prod, mykey->decr_prod, 
		      gbp->gs[(n-main_index)+receiver]);
	}
	if(0 && DEBUG) 
	  printf("added index = %d\n", main_index);
      }
      main_index++;
      working = working>>1;
    }
  }
}


  
void Change_encr_prod_indicies(global_broadcast_params_t gbp, 
			       broadcast_system_t sys, int *adds, 
			       int N_adds, int *rems, int N_rems)
{
  // REMOVES THE OLD ONES, THEN ADDS THE NEW ONES
  int i;
  element_t temp_inv;
  int incl_num;
  if(!gbp) {
    printf("ACK!  You gave me no broadcast params!  I die.\n");
    return;
  }
  if(!sys) {
    printf("ACK!  You gave me no broadcast system!  I die.\n");
    return;
  }
  int n = gbp->num_users;
  element_init(temp_inv, gbp->pairing->G1);
  if(rems) {
    for(i = 0; i < N_rems; i++) {
      //removing elements from the set after checking if it's in the bit-vector
      incl_num = rems[i];
      if(incl_num < 1 || incl_num > gbp->num_users) {
	printf("element %d was outside the range of valid users\n",i);
	printf("only give me valid values.  i die.\n");
	return;
      }
      element_invert(temp_inv, gbp->gs[n-incl_num]);
      element_mul(sys->encr_prod, sys->encr_prod, temp_inv);
    } 
  }
  if(adds) {
    for(i = 0; i < N_adds; i++) {
      //adding elements to the set after checking if they're already there
      incl_num = adds[i];
      if(incl_num < 1 || incl_num > gbp->num_users) {
	printf("element %d was outside the range of valid users\n",i);
	printf("only give me valid values.  i die.\n");
	return;
      }
      element_mul(sys->encr_prod, sys->encr_prod, gbp->gs[n-incl_num]);      
    }  
  }
}

void PrintBitString(char *bs, int length) 
{
  if(!bs) {
    printf("the bitstring you provided was null.\n");
    return;
  }
  int i,j;
  char working;
  printf("\nbitstring = ");
  for(i = 0; i < length; i++) {
    working = bs[i];
    for(j = 0; j < 8; j++) {
      if(working & 1) {
	printf("1");
      } else {
	printf("0");
      }
      working = working >> 1;
    }
  }
  printf("\n");
}


void Gen_encr_prod_from_indicies(global_broadcast_params_t gbp, 
				 broadcast_system_t sys,
				 int *in_recip, int num_recip)
{
  if(!gbp) {
    printf("ACK!  You gave me no broadcast params!  I die.\n");
    return;
  }
  if(!sys) {
    printf("ACK!  You gave me no broadcast system!  I die.\n");
    return;
  }
  if(!in_recip) {
    printf("ACK!  You gave me no recipient list!  I die.\n");
    return;
  }
  //RECIPIENT NUMBERS MUST BE IN NORMAL NOTATION e.g. 1 to N
  int i;
  element_init(sys->encr_prod, gbp->pairing->G1);
  int n = gbp->num_users;

  //UPDATE THE BIT VECTOR
  int incl_num = in_recip[0];
  if(incl_num < 1 || incl_num > n) {
    printf("element was outside the range of valid users\n");
    printf("only give me valid values.  i die.\n");
    return;
  }
  if(DEBUG && 0) printf("\nputting element %d in product",incl_num);
  element_set(sys->encr_prod, gbp->gs[n-incl_num]);
  for(i = 1; i < num_recip; i++) {
    incl_num = in_recip[i];
    if(incl_num < 1 || incl_num > n) {
      printf("element %d was outside the range of valid users\n",i);
      printf("only give me valid values.  i die.\n");
      return;
    }
    if(DEBUG && 0) printf("\nputting element %d in product",incl_num);
    element_mul(sys->encr_prod, sys->encr_prod, gbp->gs[n-incl_num]);
  }
}


void Gen_encr_prod_from_bitvec(global_broadcast_params_t gbp, 
			       broadcast_system_t sys, char *recip)
{
  // recip must have length num_users/8;
  // BITVector[0    ]lsb corresponds to Recipient 1, and
  // BITVector[n/8-1]msb corresponds to Recipient n
  char working;
  int main_index = 1;
  int i,j;
  int already_set = 0;

  if(!gbp) {
    printf("ACK!  You gave me no broadcast params!  I die.\n");
    return;
  }
  if(!sys) {
    printf("ACK!  You gave me no broadcast system!  I die.\n");
    return;
  }
  if(!recip) {
    printf("ACK!  You gave me no recipient list!  I die.\n");
    return;
  }
  element_init(sys->encr_prod, gbp->pairing->G1);

  int n = gbp->num_users;
  int length = n / 8;

  for(i = 0; i < length; i++) {
    working = recip[i];
    for(j = 0; j < 8; j++) {
      if(working & 1) {
	if(!already_set) {
	  element_set(sys->encr_prod, gbp->gs[n-main_index]);
	  already_set = 1;
	} else {
	  element_mul(sys->encr_prod, sys->encr_prod, gbp->gs[n-main_index]);
	}
      }
      main_index++;
      working = working>>1;
    }
  }
}


void Get_priv_key(global_broadcast_params_t gbp, broadcast_system_t sys,
		      int i, priv_key_t mykey)
{
  if(!gbp) {
    printf("ACK!  You gave me no broadcast params!  I die.\n");
    return;
  }
  if(!sys) {
    printf("ACK!  You gave me no broadcast system!  I die.\n");
    return;
  }
  if(!mykey) {
    printf("ACK!  You gave me no private key!  I die.\n");
    return;
  }
  if(i < 1 || i > gbp->num_users) {
    printf("ACK!  You gave me an index that's out of bounds!  I die.\n");
    printf("You must you standard notation [1...n]\n");
    return;    
  }
  element_init(mykey->g_i_gamma, gbp->pairing->G1);
  element_init(mykey->g_i, gbp->pairing->G1);
  element_init(mykey->h_i, gbp->pairing->G2);
  mykey->index = i;
  element_set(mykey->g_i,gbp->gs[i-1]);
  element_set(mykey->h_i,gbp->hs[i-1]);
  element_pow(mykey->g_i_gamma, gbp->gs[i-1],sys->priv_key);
}



void Gen_broadcast_system(global_broadcast_params_t gbp,
			  broadcast_system_t *sys)
{
  if(!gbp) {
    printf("ACK!  You gave me no broadcast params!  I die.\n");
    return;
  }
  broadcast_system_t my_sys;
  my_sys = malloc(sizeof(struct broadcast_system_s));

  mpz_init(my_sys->priv_key);
  
  pbc_mpz_random(my_sys->priv_key, gbp->pairing->r);
  
  if(mpz_sgn(my_sys->priv_key) < 0) {
    mpz_neg(my_sys->priv_key, my_sys->priv_key);
  }
  element_init(my_sys->pub_key, gbp->pairing->G1);
  element_pow(my_sys->pub_key, gbp->g, my_sys->priv_key); 

  *sys = my_sys;
}


void Setup_global_broadcast_params(global_broadcast_params_t *sys, 
				   int num_users, char *pairFileName)
{
  global_broadcast_params_t gbs;

  gbs = malloc(sizeof(struct global_broadcast_params_s));
  
  // Setup curve in gbp

  FILE *curveFile = fopen(pairFileName, "r");
  gbs->pairFileName = strdup(pairFileName);
  if(!curveFile) {
    printf("%s doesn't exist!  exiting! \n\n", pairFileName);
    return;
  }
  
  pairing_init_inp_str(gbs->pairing, curveFile);
  fclose(curveFile);

  gbs->num_users = num_users;
  element_t *lgs;
  element_t *lhs;
  int i;

  if(num_users % 8 != 0) {
    printf("\nSystem size must be a multiple of 8\n");
    printf("Didn't finish system setup\n\n");
    return;
  }
  lgs = malloc(2 * num_users * sizeof(element_t));
  lhs = malloc(2 * num_users * sizeof(element_t));
  if(!(lhs) || !(lgs)) {
    printf("\nMalloc Failed\n");
    printf("Didn't finish system setup\n\n");
  }
  //Choosing random G & H
  element_init(gbs->g, gbs->pairing->G1);
  element_random(gbs->g);
  element_init(gbs->h, gbs->pairing->G2);
  element_random(gbs->h);

  mpz_t alpha;
  //Pick a random exponent alpha
  mpz_init(alpha);
  pbc_mpz_random(alpha, gbs->pairing->r);
  if(mpz_sgn(alpha) < 0) {
    mpz_neg(alpha,alpha);
  }

  //Make the 0th elements equal to x^alpha
  element_init(lgs[0], gbs->pairing->G1);
  element_init(lhs[0], gbs->pairing->G2);
  element_pow(lgs[0],gbs->g, alpha);
  element_pow(lhs[0],gbs->h, alpha);

  //Fill in the gs and the hs arrays
  for(i = 1; i < 2*num_users; i++) { 
    //raise alpha to one more power
    if(DEBUG) {
      if(!(i % 5)) 
	printf("Finished computing elem %d\n",i);
    }

    element_init(lgs[i], gbs->pairing->G1);
    element_pow(lgs[i],lgs[i-1], alpha);
    element_init(lhs[i], gbs->pairing->G2);
    element_pow(lhs[i],lhs[i-1], alpha);
    if(i == num_users+1) {
      element_clear(lgs[i-1]);
      element_clear(lhs[i-1]);
    }
	   
  }
  
  //For simplicity & so code was easy to read
  gbs->gs = lgs;
  gbs->hs = lhs;
  
  *sys = gbs;
  mpz_clear(alpha);
}



