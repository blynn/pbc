#include "pbc.h"
#include "pbc_test.h"

int main(void) {
  element_t g, h;
  element_t x1, x2;
  element_t zg, zh, z;
  pairing_t pairing;
  pbc_param_t param;

  pbc_param_init_a_gen(param, 160, 512);
  pairing_init_pbc_param(pairing, param);
  pbc_param_clear(param);

  element_init_G1(g, pairing);
  element_init_G1(zg, pairing);
  element_init_G2(h, pairing);
  element_init_G2(zh, pairing);
  element_init_GT(x1, pairing);
  element_init_GT(x2, pairing);
  element_init_Zr(z, pairing);
  element_random(g);
  element_random(h);
  element_printf("g = %B\n", g);
  element_printf("h = %B\n", h);
  pairing_apply(x1, g, h, pairing);
  element_printf("f(g, h) = %B\n", x1);

  element_random(z);
  element_printf("z = %B\n", z);

  element_pow_zn(x1, x1, z);
  element_printf("f(g, h)^z = %B\n", x1);

  element_pow_zn(zg, g, z);
  element_printf("g^z = %B\n", zg);
  pairing_apply(x2, zg, h, pairing);
  element_printf("f(g^z, h) = %B\n", x2);

  EXPECT(!element_cmp(x1, x2));

  element_pow_zn(zh, h, z);
  element_printf("h^z = %B\n", zh);
  pairing_apply(x2, g, zh, pairing);
  element_printf("f(g, h^z) = %B\n", x2);

  EXPECT(!element_cmp(x1, x2));

  {
    int i;
    int len = element_length_in_bytes(h);
    printf("length_in_bytes(h) = %d\n", len);
    unsigned char *data = pbc_malloc(len);
    element_to_bytes(data, h);
    for (i=0; i<len; i++) {
      printf(" %02X", data[i]);
      if (15 == (i % 16)) printf("\n");
    }
    printf("\n");
    element_from_bytes(h, data);
    element_printf("from_bytes h = %B\n", h);
  }

  element_clear(g);
  element_clear(h);
  element_clear(x1);
  element_clear(x2);
  element_clear(zg);
  element_clear(zh);
  element_clear(z);
  pairing_clear(pairing);
  return 0;
}
