// Boneh-Lynn-Shacham short signatures demo.
//
// See the PBC_sig library for a practical implementation.
//
// Ben Lynn
#include <pbc.h>
#include <pbc_test.h>

int main(int argc, char **argv) {
  pairing_t pairing;
  element_t g, h;
  element_t public_key, sig;
  element_t secret_key;
  element_t temp1, temp2;

  pbc_demo_pairing_init(pairing, argc, argv);

  element_init_G2(g, pairing);
  element_init_G2(public_key, pairing);
  element_init_G1(h, pairing);
  element_init_G1(sig, pairing);
  element_init_GT(temp1, pairing);
  element_init_GT(temp2, pairing);
  element_init_Zr(secret_key, pairing);

  printf("Short signature test\n");

  //generate system parameters
  element_random(g);
  element_printf("system parameter g = %B\n", g);

  //generate private key
  element_random(secret_key);
  element_printf("private key = %B\n", secret_key);

  //compute corresponding public key
  element_pow_zn(public_key, g, secret_key);
  element_printf("public key = %B\n", public_key);

  //generate element from a hash
  //for toy pairings, should check that pairing(g, h) != 1
  element_from_hash(h, "hashofmessage", 13);
  element_printf("message hash = %B\n", h);

  //h^secret_key is the signature
  //in real life: only output the first coordinate
  element_pow_zn(sig, h, secret_key);
  element_printf("signature = %B\n", sig);

  {
    int n = pairing_length_in_bytes_compressed_G1(pairing);
    //int n = element_length_in_bytes_compressed(sig);
    int i;
    unsigned char *data = pbc_malloc(n);

    element_to_bytes_compressed(data, sig);
    printf("compressed = ");
    for (i = 0; i < n; i++) {
      printf("%02X", data[i]);
    }
    printf("\n");

    element_from_bytes_compressed(sig, data);
    element_printf("decompressed = %B\n", sig);

    pbc_free(data);
  }

  //verification part 1
  element_pairing(temp1, sig, g);
  element_printf("f(sig, g) = %B\n", temp1);

  //verification part 2
  //should match above
  element_pairing(temp2, h, public_key);
  element_printf("f(message hash, public_key) = %B\n", temp2);

  if (!element_cmp(temp1, temp2)) {
    printf("signature verifies\n");
  } else {
    printf("*BUG* signature does not verify *BUG*\n");
  }

  {
    int n = pairing_length_in_bytes_x_only_G1(pairing);
    //int n = element_length_in_bytes_x_only(sig);
    int i;
    unsigned char *data = pbc_malloc(n);

    element_to_bytes_x_only(data, sig);
    printf("x-coord = ");
    for (i = 0; i < n; i++) {
      printf("%02X", data[i]);
    }
    printf("\n");

    element_from_bytes_x_only(sig, data);
    element_printf("de-x-ed = %B\n", sig);

    element_pairing(temp1, sig, g);
    if (!element_cmp(temp1, temp2)) {
      printf("signature verifies on first guess\n");
    } else {
      element_invert(temp1, temp1);
      if (!element_cmp(temp1, temp2)) {
        printf("signature verifies on second guess\n");
      } else {
        printf("*BUG* signature does not verify *BUG*\n");
      }
    }

    pbc_free(data);
  }

  //a random signature shouldn't verify
  element_random(sig);
  element_pairing(temp1, sig, g);
  if (element_cmp(temp1, temp2)) {
    printf("random signature doesn't verify\n");
  } else {
    printf("*BUG* random signature verifies *BUG*\n");
  }

  element_clear(sig);
  element_clear(public_key);
  element_clear(secret_key);
  element_clear(g);
  element_clear(h);
  element_clear(temp1);
  element_clear(temp2);
  pairing_clear(pairing);
  return 0;
}
