/* test ternary extension fields $GF(3^m)$, $GF(3^{2*m})$, $GF(3^{3*m})$ and $GF(3^{6*m})$
   Outputing nothing if everything is good. */

#include "pbc.h"
#include "pbc_ternary_extension_field.h"
#include "pbc_test.h"
#include <string.h>
#include <stdio.h>

typedef struct {
    unsigned int len;
    unsigned int m;
    unsigned int t;
    element_ptr p;
} params;

#define data(x) ((unsigned long*)x->data)
#define params(x) ((params *)x->field->data)
#define print(e) {printf(#e": "); element_out_str(stdout, 0, e); printf("\n");}

static field_t f97, f97_2, f97_3, f97_6;
static element_t e0, e1, e2, a, b, a2, b2, a3, b3, a6, b6;
static unsigned char *data;

static void test_gf3m_param(void) {
    params *pa = (params *) f97->data;
    element_to_bytes(data, pa->p);
    unsigned i;
    unsigned char w;
    for (i = 0; i < pa->len * 2 * sizeof(unsigned long); i++) {
        switch (i) {
        case 1:
            w = 1;
            break; // 2
        case 2:
            w = 16;
            break; // x^12
        case 24:
            w = 2;
            break; // x^97
        default:
            w = 0;
        }
        EXPECT(data[i] == w);
    }
}

static void test_gf3m_to_bytes(void) {
    element_random(a);
    element_to_bytes(data, a);
    element_from_bytes(b, data);
    EXPECT(0 == element_cmp(a, b));
}

static void test_gf3m_add(void) {
    element_random(a);
    element_add(b, a, a);
    element_add(b, b, b);
    element_sub(b, b, a);
    element_sub(b, b, a);
    element_sub(b, b, a);
    EXPECT(!element_cmp(a, b));

    element_add(b, params(a)->p, a);
    element_sub(b, b, params(a)->p);
    EXPECT(!element_cmp(a, b));
}

static void test_gf3m_neg(void) {
    element_random(a);
    element_neg(b, a);
    element_add(b, a, b);
    EXPECT(!element_cmp(b, e0));
}

static void test_gf3m_mult(void) {
    element_random(a);
    element_mul(a, a, e0);
    EXPECT(!element_cmp(a, e0));

    element_random(a);
    element_mul(b, a, e1);
    EXPECT(!element_cmp(a, b));

    element_random(a);
    element_mul(b, a, e2);
    element_add(a, a, b);
    EXPECT(!element_cmp(a, e0));
}

static void test_gf3m_cubic(void) {
    element_random(a);
    element_mul(b, a, a);
    element_mul(b, a, b);
    element_cubic(a, a);
    EXPECT(!element_cmp(a, b));
}

static void test_gf3m_cubic2(void) {
    unsigned long x[] = {1153286547535200267ul, 6715371622ul, 4990694927524257316ul,  210763913ul};
    unsigned long y[] = {8145587063258678275ul, 6451025920ul, 9976895054123379152ul, 1275593166ul};
    memcpy(a->data, x, sizeof(x));
    memcpy(b->data, y, sizeof(y));
    element_cubic(a, a);
    EXPECT(!element_cmp(a, b));
}

static void test_gf3m_inverse(void) {
    element_set1(a);
    element_invert(b, a);
    EXPECT(!element_cmp(b, e1));

    element_set(a, e2);
    element_invert(b, a);
    EXPECT(!element_cmp(b, e2));

    element_random(a);
    element_invert(b, a);
    element_mul(a, a, b);
    EXPECT(!element_cmp(a, e1));
}

static void test_gf3m_sqrt(void) {
    mpz_t t;
    mpz_init(t);
    mpz_sub_ui(t, a->field->order, 1); // t == field_order - 1
    element_random(a);
    element_pow_mpz(a, a, t);
    EXPECT(!element_cmp(a, e1));

    while(1){
        element_random(a);
        element_mul(b, a, a);
        element_sqrt(b, b);
        if(element_cmp(a, b)) {// a != b
            element_neg(b, b);
            if(!element_cmp(a, b)) break;
        }
    }
    mpz_clear(t);
}

static void test_gf32m_cubic(void) {
    element_random(a2);
    element_mul(b2, a2, a2);
    element_mul(b2, b2, a2);
    element_cubic(a2, a2);
    EXPECT(!element_cmp(a2, b2));
}

static void test_gf33m_cubic(void) {
    element_random(a3);
    element_mul(b3, a3, a3);
    element_mul(b3, b3, a3);
    element_cubic(a3, a3);
    EXPECT(!element_cmp(a3, b3));
}

static void test_gf33m_inverse(void) {
    element_random(a3);
    element_invert(b3, a3);
    element_mul(a3, a3, b3);
    element_ptr a0 = element_item(a3, 0);
    EXPECT(!element_cmp(a0, e1));
}

static void test_gf36m_cubic(void) {
    element_random(a6);
    element_mul(b6, a6, a6);
    element_mul(b6, b6, a6);
    element_cubic(a6, a6);
    EXPECT(!element_cmp(a6, b6));
}

void setup(void) {
    field_init_gf3m(f97, 97, 12);
    element_init(a, f97);
    element_init(b, f97);
    element_init(e0, f97);
    element_init(e1, f97);
    element_init(e2, f97);
    element_set1(e1);
    element_neg(e2, e1);

    field_init_gf32m(f97_2, f97);
    element_init(a2, f97_2);
    element_init(b2, f97_2);

    field_init_gf33m(f97_3, f97);
    element_init(a3, f97_3);
    element_init(b3, f97_3);

    field_init_gf33m(f97_6, f97_2);
    element_init(a6, f97_6);
    element_init(b6, f97_6);

    data = pbc_malloc(f97->fixed_length_in_bytes);
}

void tear_down(void) {
    pbc_free(data);

    element_clear(e0);
    element_clear(e1);
    element_clear(e2);
    element_clear(a);
    element_clear(b);
    element_clear(a2);
    element_clear(b2);
    element_clear(a3);
    element_clear(b3);
    element_clear(a6);
    element_clear(b6);

    field_clear(f97_6);
    field_clear(f97_3);
    field_clear(f97_2);
    field_clear(f97);
}

int main(void) {
    setup();

    test_gf3m_param();
    test_gf3m_to_bytes();
    test_gf3m_add();
    test_gf3m_neg();
    test_gf3m_mult();
    test_gf3m_cubic();
    test_gf3m_cubic2();
    test_gf3m_inverse();
    test_gf3m_sqrt();
    test_gf32m_cubic();
    test_gf33m_cubic();
    test_gf33m_inverse();
    test_gf36m_cubic();

    tear_down();
    return 0;
}
