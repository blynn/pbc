//Boneh-Franklin Identity-Based Encryption demo
//Note: in real life it may be better to swap the roles of G1 and G2
//although the system parameters take longer to compute
//(and more room to store), hashing ID's to elements of G1
//is faster than to G2, and this has to be done often.
#include "pbc.h"

int main(void)
{
    element_t g, h, s;
    element_t rg, zg, zh;
    pairing_t pairing;
    mpz_t master, r;


    pairing_init_inp_str(pairing, stdin);
    element_init_G1(g, pairing);
    element_init_G1(zg, pairing);
    element_init_G1(rg, pairing);
    element_init_G2(h, pairing);
    element_init_G2(zh, pairing);
    element_init_GT(s, pairing);
    mpz_init(master);
    mpz_init(r);

    printf("Identity-based encryption test program\n");

    //generate master secret
    pbc_mpz_random(master, pairing->r);
    printf("master secret = ");
    mpz_out_str(stdout, 0, master);
    printf("\n");

    //generate g, compute g^master
    element_random(g);

    printf("g = ");
    element_out_str(stdout, 0, g);
    printf("\n");
    element_pow(zg, g, master);
    printf("g^master = ");
    element_out_str(stdout, 0, zg);
    printf("\n");

    //pick random h, which represents what an ID might hash to
    //for toy examples, should check that pairing(g, h) != 1
    element_random(h);
    printf("ID hashes to ");
    element_out_str(stdout, 0, h);
    printf("\n");

    //h^master is the corresponding private key
    element_pow(zh, h, master);
    printf("private key = ");
    element_out_str(stdout, 0, zh);
    printf("\n");

    //encryption: first pick random r
    pbc_mpz_random(r, pairing->r);
    printf("random r = ");
    mpz_out_str(stdout, 0, r);
    printf("\n");

    //compute s = f(g^master, h)^r, used to encrypt the message
    bilinear_map(s, zg, h, pairing);
    element_pow(s, s, r);
    printf("f(g^master, h)^r = ");
    element_out_str(stdout, 0, s);
    printf("\n");

    //we transmit g^r along with the encryption
    element_pow(rg, g, r);
    printf("g^r = ");
    element_out_str(stdout, 0, rg);
    printf("\n");

    //decryption: compute f(g^r, h^master)
    //should equal s
    bilinear_map(s, rg, zh, pairing);
    printf("f(g^r, h^master) = ");
    element_out_str(stdout, 0, s);
    printf("\n");

    return 0;
}
