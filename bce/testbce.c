/* Implementation of Boneh-Waters broadcast encryption scheme
   Code by:  Matt Steiner   MattS@cs.stanford.edu
   testbce.c
*/

#include "pbc.h"
#include <gmp.h>
#include <string.h>
#include "bce.h"

#define N 64
#define N_DIV_EIGHT  N/8


int main(void)
{
  int i;
  
  global_broadcast_params_t gbs;

  //Global Setup
  Setup_global_broadcast_params(&gbs, N, "c201.param");
  
  if(0 && DEBUG) {
    printf("\ng = ");  
    element_out_str(stdout, 0, gbs->g);
    printf("\nh = ");  
    element_out_str(stdout, 0, gbs->h);
    for(i = 0; i < 1; i++) {
      printf("\nThe next element is %d------------------------------------",i);
      printf("\ngs[%d] = ", i);
      element_out_str(stdout, 0, gbs->gs[i]);
      printf("\nhs[%d] = ",i);
      element_out_str(stdout, 0, gbs->hs[i]);
    }
    printf("\n");
  }
  
  //Broadcast System Setup
  broadcast_system_t sys;
  Gen_broadcast_system(gbs, &sys);
  
  struct single_priv_key_s mykey;
  struct single_priv_key_s mykey2;
  struct single_priv_key_s mykey3;
  
  Get_priv_key(gbs, sys, 2, &mykey);
  //if(DEBUG) printf("done 1\n");
  //if(DEBUG) printf("done 2\n");
  Get_priv_key(gbs, sys, 2, &mykey3);
  //if(DEBUG) printf("done 3\n");
  
  if(DEBUG && 0) {
    printf("\ng_i = ");
    element_out_str(stdout, 0, mykey.g_i);
    printf("\nh_i = ");
    element_out_str(stdout, 0, mykey.h_i);
    printf("\ng_i_gamma = ");
    element_out_str(stdout, 0, mykey.g_i_gamma);
    printf("\n");
    printf("\ng_i = ");
    element_out_str(stdout, 0, mykey2.g_i);
    printf("\nh_i = ");
    element_out_str(stdout, 0, mykey2.h_i);
    printf("\ng_i_gamma = ");
    element_out_str(stdout, 0, mykey2.g_i_gamma);
    printf("\n");
     printf("\ng_i = ");
    element_out_str(stdout, 0, mykey3.g_i);
    printf("\nh_i = ");
    element_out_str(stdout, 0, mykey3.h_i);
    printf("\ng_i_gamma = ");
    element_out_str(stdout, 0, mykey3.g_i_gamma);
    printf("\n");
 }

  char recip[N_DIV_EIGHT];
  for(i = 0; i < 2; i++) recip[i] = 254;
  for(i = 2; i < N_DIV_EIGHT; i++) recip[i] = 0;

  Gen_encr_prod_from_bitvec(gbs, sys, recip);
  //Product_Is_Right(gbs, sys, recip);
  //TESTING FOR SYSTEM LOAD AND STORE
  global_broadcast_params_t gbp2;
  broadcast_system_t sys2;
  global_broadcast_params_t gbp3;
  broadcast_system_t sys3;

  StoreParams("system.stor", gbs, sys);
  //printf("\ndone storing!!!!!!!!!\n\n");
  LoadParams("system.stor", &gbp2, &sys2);
  LoadParams("system.stor", &gbp3, &sys3);

  //printf("\ndone loading!!!!!!!!!\n\n");
  //StoreParams("system2.stor", "pairing2.stor", gbp2, sys2);
  //LoadParams("system2.stor", "pairing2.stor", &gbs, &sys);

  Get_priv_key(gbs, sys, 2, &mykey2);

  if(DEBUG) {
    printf("\noldg = ");
    element_out_str(stdout, 0, gbs->g);
    printf("\nnew = ");
    element_out_str(stdout, 0, gbp2->g);
    printf("\noldh = ");
    element_out_str(stdout, 0, gbs->h);
    printf("\nnew = ");
    element_out_str(stdout, 0, gbp2->h);
    printf("\noldgs = ");
    element_out_str(stdout, 0, gbs->gs[0]);
    printf("\nnew = ");
    element_out_str(stdout, 0, gbp2->gs[0]);
    printf("\nold = ");
    element_out_str(stdout, 0, gbs->gs[31]);
    printf("\nnew = ");
    element_out_str(stdout, 0, gbp2->gs[31]);
    printf("\noldhs = ");
    element_out_str(stdout, 0, gbs->hs[0]);
    printf("\nnew = ");
    element_out_str(stdout, 0, gbp2->hs[0]);
    printf("\nold = ");
    element_out_str(stdout, 0, gbs->hs[31]);
    printf("\nnew = ");
    element_out_str(stdout, 0, gbp2->hs[31]);
    printf("\n old n_u = %d", gbs->num_users);
    printf("\n new n_u = %d", gbp2->num_users); 
    printf("\nolde = ");
    element_out_str(stdout, 0, sys->encr_prod);
    printf("\nnew = ");
    element_out_str(stdout, 0, sys2->encr_prod);
    printf("\noldp = ");
    element_out_str(stdout, 0, sys->pub_key);
    printf("\nnew = ");
    element_out_str(stdout, 0, sys2->pub_key);
  }

  
  //int in_recip[5] = {4, 5, 6, 7, 8 };
  //int num_recip = 5;
  //int rems[3] = { 5, 6, 7 };
  //int N_rems = 3;
  //int adds[12] = { 2, 3, 5, 6, 7, 10, 11, 12, 13, 14, 15, 16 };
  //int N_adds = 12;
  // FINAL ELEMENTS IN PRODUCT SHOULD BE 2-8, & 10-16

  /*
  Gen_encr_prod_from_indicies(gbs, sys2, in_recip, num_recip);

  if(DEBUG) {
    PrintBitString(sys2->recipients,BSL);
    printf("\nsys2 encr_product = ");
    element_out_str(stdout, 0, sys2->encr_prod);
    printf("\n");
  }

  Change_encr_prod_indicies(gbs, sys2, adds, N_adds, rems, N_rems);
  if(DEBUG) {
    PrintBitString(sys2->recipients,BSL);
    printf("\nsys2 encr_product = ");
    element_out_str(stdout, 0, sys2->encr_prod);
    printf("\n");
  }
    

  if(DEBUG) {
    PrintBitString(sys->recipients,BSL);
    printf("\nsys1 encr_product = ");
    element_out_str(stdout, 0, sys->encr_prod);
  }  
  */
  
  Gen_decr_prod_from_bitvec(gbs, 2, recip, &mykey);
  //if(DEBUG && 0) printf("\ndone 1 decr\n");
  Gen_decr_prod_from_bitvec(gbs, 2, recip, &mykey2);
  //if(DEBUG && 0) printf("\ndone 2 decr\n");
  Gen_decr_prod_from_bitvec(gbs, 2, recip, &mykey3);
  //if(DEBUG && 0) printf("\ndone 3 decr\n");
  //Gen_decr_prod_from_indicies(gbs, 2, in_recip, num_recip, &mykey2);  
  //Change_decr_prod_indicies(gbs, 2, adds, N_adds, rems, N_rems, &mykey2);

  //Gen_decr_prod_from_bitvec(gbs, 2, recip, &mykey3);
 

  if(0 && DEBUG) {
    printf("\n");
    printf("mykey1 decr_product = ");
    element_out_str(stdout, 0, mykey.decr_prod);
    printf("\n");
  }  
  if(DEBUG && 0) {
    printf("\n");
    printf("mykey2 decr_product = ");
    element_out_str(stdout, 0, mykey2.decr_prod);
    printf("\n");
  }
  if(DEBUG && 0) {
    printf("\n");
    printf("mykey3 decr_product = ");
    element_out_str(stdout, 0, mykey3.decr_prod);
    printf("\n");
  }
 



  //TESTING FOR SINGLE KEY LOAD AND STORE
  priv_key_t load_key = (priv_key_t)malloc(sizeof(struct single_priv_key_s));

  StorePrivKey("key2.stor", &mykey);
  LoadPrivKey("key2.stor", &load_key, gbs);
  
  if(DEBUG) {
    printf("\nold = ");
    element_out_str(stdout, 0, mykey.g_i_gamma);
    printf("\nnew = ");
    element_out_str(stdout, 0, load_key->g_i_gamma);
    printf("\nold = ");
    element_out_str(stdout, 0, mykey.g_i);
    printf("\nnew = ");
    element_out_str(stdout, 0, load_key->g_i);
    printf("\nold = ");
    element_out_str(stdout, 0, mykey.h_i);
    printf("\nnew = ");
    element_out_str(stdout, 0, load_key->h_i);
    printf("\nold = ");
    element_out_str(stdout, 0, mykey.decr_prod);
    printf("\nnew = ");
    element_out_str(stdout, 0, load_key->decr_prod);
    printf("\n index = %d", mykey.index);
    printf("\n index = %d", load_key->index); 
  }

  ct_t myCT = (ct_t) malloc(sizeof(struct ciphertext_s));
  ct_t myCT2 = (ct_t) malloc(sizeof(struct ciphertext_s));
  ct_t myCT3 = (ct_t) malloc(sizeof(struct ciphertext_s));
  //int recip2[14] = { 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16 };
  //int n_recip2 = 14; 
  element_t key1;
  element_t key2;
  element_t key3;
  element_t key4;
  element_t key5;
  element_t key6;

  BroadcastKEM_using_product(gbs, sys, myCT, key1);
  DecryptKEM_using_product(gbs, &mykey, key4, myCT);
  BroadcastKEM_using_product(gbs, sys, myCT3, key3);
  DecryptKEM_using_product(gbp3, &mykey3, key6, myCT3);
  BroadcastKEM_using_product(gbs, sys, myCT2, key2);
  DecryptKEM_using_product(gbp2, &mykey2, key5, myCT2);


  //BroadcastKEM_using_bitvec(gbs, sys, recip, myCT2, key2);
  //BroadcastKEM_using_indicies(gbs, sys, myCT3, recip2, n_recip2, key3);

  
  if(DEBUG) {
    //COMPARE ALL THREE CTs!
    printf("\n1-C0 = ");
    element_out_str(stdout, 0, myCT->C0);
    printf("\n2-C0 = ");
    element_out_str(stdout, 0, myCT2->C0);
    printf("\n3-C0 = ");
    element_out_str(stdout, 0, myCT3->C0);
    printf("\n1-C1 = ");
    element_out_str(stdout, 0, myCT->C1);
    printf("\n2-C1 = ");
    element_out_str(stdout, 0, myCT2->C1);
    printf("\n3-C1 = ");
    element_out_str(stdout, 0, myCT3->C1);
  }
  
  
  printf("\nkey1 = ");
  element_out_str(stdout, 0, key1);
  printf("\n");
  printf("\nkey2 = ");
  element_out_str(stdout, 0, key2);
  printf("\n");
  printf("\nkey3 = ");
  element_out_str(stdout, 0, key3);
  printf("\n");

  //PrintBitString(mykey.recipients, BSL);
  //DecryptKEM_using_product(gbs, &mykey2, key5, myCT2);
  

  //printf("\nmyprivkey = ");
  //element_out_str(stdout, 0, mykey.g_i_gamma);
  //printf("\n");
  printf("\nkey1 = ");
  element_out_str(stdout, 0, key4);
  printf("\n");
  printf("\nkey2 = ");
  element_out_str(stdout, 0, key5);
  printf("\n");
  printf("\nkey3 = ");
  element_out_str(stdout, 0, key6);
  printf("\n");

  FreeCT(myCT);
  FreeBCS(sys);
  FreeGBP(gbs);
  FreeGBP(gbp2);
  FreeBCS(sys2);
  FreePK(&mykey);
  return 0;

}




