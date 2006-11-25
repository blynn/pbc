// Boneh-Franklin Identity-Based Encryption demo
//
// Note: in real life it may be better to swap the roles
// of G1 and G2 for some pairing types. Although the system
// parameters take longer to compute (and more room to store),
// hashing ID's to elements of G1 is faster than to G2, and
// this has to be done often.
//
// Ben Lynn
#include "pbc.h"

int main(void)
{
    element_t g, h, s;
    element_t rg, zg, zh;
    pairing_t pairing;
    element_t master, r;

    pairing_init_inp_str(pairing, stdin);
    element_init_G1(g, pairing);
    element_init_G1(zg, pairing);
    element_init_G1(rg, pairing);
    element_init_G2(h, pairing);
    element_init_G2(zh, pairing);
    element_init_GT(s, pairing);
    element_init_Zr(master, pairing);
    element_init_Zr(r, pairing);

    printf("Identity-based encryption test program\n");

    //generate master secret
    element_random(master);
    element_printf("master secret = %B\n", master);

    //generate g, compute g^master
    element_random(g);

    element_printf("g = %B\n", g);
    element_pow_zn(zg, g, master);
    element_printf("g^master = %B\n", zg);

    //pick random h, which represents what an ID might hash to
    //for toy examples, should check that pairing(g, h) != 1
    element_random(h);
    element_printf("ID hashes to %B\n", h);

    //h^master is the corresponding private key
    element_pow_zn(zh, h, master);
    element_printf("private key = %B\n", zh);

    //encryption: first pick random r
    element_random(r);
    element_printf("random r = %B\n", r);

    //compute s = f(g^master, h)^r, used to encrypt the message
    pairing_apply(s, zg, h, pairing);
    element_pow_zn(s, s, r);
    element_printf("f(g^master, h)^r = %B\n", s);

    //we transmit g^r along with the encryption
    element_pow_zn(rg, g, r);
    element_printf("g^r = %B\n", rg);

    //decryption: compute f(g^r, h^master)
    //should equal s
    pairing_apply(s, rg, zh, pairing);
    element_printf("f(g^r, h^master) = %B\n", s);

    element_clear(g);
    element_clear(h);
    element_clear(s);
    element_clear(rg);
    element_clear(zg);
    element_clear(zh);
    element_clear(master);
    element_clear(r);
    pairing_clear(pairing);
    return 0;
}
