/* Test eta_T pairing over ternary extension fields.
   Outputing nothing if everything is good. */

#include <stddef.h>
#include <stdarg.h>
#include <stdio.h>
#include <gmp.h>
#include "pbc.h"
#include "pbc_ternary_extension_field.h"
#include "pbc_test.h"

static pairing_t pairing;
static element_t a1, a2, b1, b2, c1, c2;
static mpz_t order;

static void setup(void) {
    mpz_init(order);
    mpz_set_str(order, "2726865189058261010774960798134976187171462721", 10);
    const char *param = "type i\n" "m 97\n" "t 12\n" "n2 7\n"
                        "n 2726865189058261010774960798134976187171462721\n";
    EXPECT(pairing_init_set_str(pairing, param) == 0);
    element_init_G1(a1, pairing);
    element_init_G1(a2, pairing);
    element_init_G2(b1, pairing);
    element_init_G2(b2, pairing);
    element_init_GT(c1, pairing);
    element_init_GT(c2, pairing);
}

static void test_set_mpz(void) {
    mpz_t a;
    mpz_init(a);
    int i;
    for(i = 0; i < 2; i ++) {
        mpz_set_si(a, i);
        element_set_mpz(a1, a);
        EXPECT(element_is0(a1) && element_is1(a1));
        element_set_mpz(b1, a);
        EXPECT(element_is0(b1) && element_is1(b1));
        element_set_mpz(c1, a);
        EXPECT(element_is0(c1) && element_is1(c1));
    }
    mpz_clear(a);
}

static void test_order(void) {
    EXPECT(mpz_cmp(pairing->G1->order, order) == 0);
    EXPECT(mpz_cmp(pairing->G2->order, order) == 0);
    EXPECT(mpz_cmp(pairing->GT->order, order) == 0);
    int i;
    for (i=0; i<10; i++) {
        element_random(a1);
        EXPECT(element_is0(a1) == 0);
        element_pow_mpz(a1, a1, order);
        EXPECT(element_is0(a1));
        element_random(c1);
        EXPECT(element_is0(c1) == 0);
        element_pow_mpz(c1, c1, order);
        EXPECT(element_is0(c1));
    }
}

static void test_bilinear_with_zero(void) {
    element_set0(a1);
    element_random(b1);
    element_pairing(c1, a1, b1);
    EXPECT(element_is0(c1) && element_is1(c1));
    element_random(a1);
    element_set0(b1);
    element_pairing(c1, a1, b1);
    EXPECT(element_is0(c1) && element_is1(c1));
    element_set0(a1);
    element_set0(b1);
    element_pairing(c1, a1, b1);
    EXPECT(element_is0(c1) && element_is1(c1));
}

static void test_bilinear(void) {
    element_random(a1);
    element_mul_si(a2, a1, 33);
    element_random(b1);
    element_mul_si(b2, b1, 33);
    element_pairing(c1, a1, b2);
    element_pairing(c2, a2, b1);
    EXPECT(element_cmp(c1, c2) == 0);
    element_mul_mpz(c1, c1, order);
    EXPECT(element_is0(c1));
}

static void test_gen_param(void) {
    typedef struct {
        unsigned int len;
        int m;
        int t;
        element_ptr p;
        mpz_t n;
        mpz_t n2;
    } params;

    pbc_param_t par;
    pbc_param_init_i_gen(par, 150);
    params *p = par->data;
    EXPECT(p->m == 97);
    EXPECT(p->t == 12);
    EXPECT(!mpz_cmp(p->n, order));
    EXPECT(!mpz_cmp_ui(p->n2, 7));
    pbc_param_clear(par);
}

static void tear_down(void) {
    element_clear(a1);
    element_clear(a2);
    element_clear(b1);
    element_clear(b2);
    element_clear(c1);
    element_clear(c2);
    pairing_clear(pairing);
    mpz_clear(order);
}

int main(void) {
    setup();
    test_set_mpz();
    test_order();
    test_bilinear_with_zero();
    test_bilinear();
    test_gen_param();
    tear_down();
    return 0;
}
