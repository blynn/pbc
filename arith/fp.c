#include <ctype.h>
#include <stdio.h>
#include <gmp.h>
#include <string.h>
#include "pbc_field.h"
#include "pbc_fp.h"

void fp_tonelli(element_ptr x, element_ptr a)
{
    int s;
    int i;
    mpz_t e;
    mpz_t t, t0;
    element_t ginv, e0;
    element_ptr nqr;

    mpz_init(t);
    mpz_init(e);
    mpz_init(t0);
    element_init(ginv, a->field);
    element_init(e0, a->field);
    nqr = field_get_nqr(a->field);

    element_invert(ginv, nqr); 

    //let q be the order of the field
    //q - 1 = 2^s t, t odd
    mpz_sub_ui(t, a->field->order, 1);
    s = mpz_scan1(t, 0);
    mpz_tdiv_q_2exp(t, t, s);
    mpz_set_ui(e, 0);
    for (i=2; i<=s; i++) {
	mpz_sub_ui(t0, a->field->order, 1);
	mpz_tdiv_q_2exp(t0, t0, i);
	element_pow_mpz(e0, ginv, e);
	element_mul(e0, e0, a);
	element_pow_mpz(e0, e0, t0);
	if (!element_is1(e0)) mpz_setbit(e, i-1);
    }
    element_pow_mpz(e0, ginv, e);
    element_mul(e0, e0, a);
    mpz_add_ui(t, t, 1);
    mpz_tdiv_q_2exp(t, t, 1);
    mpz_tdiv_q_2exp(e, e, 1);

    //as suggested by Hovav Shacham, replace next three lines with
    //  element_pow2_mpz(x, e0, t, nqr, e);
    //once sliding windows for pow2 has been done
    element_pow_mpz(e0, e0, t);
    element_pow_mpz(x, nqr, e);
    element_mul(x, x, e0);

    mpz_clear(t);
    mpz_clear(e);
    mpz_clear(t0);
    element_clear(ginv);
    element_clear(e0);
}

static void (*option_fpinit)(field_ptr f, mpz_t prime) = field_init_mont_fp;

void field_init_fp(field_ptr f, mpz_t modulus)
{
    if (mpz_fits_ulong_p(modulus)) {
	field_init_naive_fp(f, modulus);
    } else {
	if (mpz_odd_p(modulus)) {
	    option_fpinit(f, modulus);
	} else {
	    field_init_faster_fp(f, modulus);
	}
    }
}

void pbc_tweak_use_fp(char *s)
{
    if (!strcmp(s, "naive")) {
	option_fpinit = field_init_naive_fp;
    } else if (!strcmp(s, "fast")) {
	option_fpinit = field_init_fast_fp;
    } else if (!strcmp(s, "faster")) {
	option_fpinit = field_init_faster_fp;
    } else if (!strcmp(s, "mont")) {
	option_fpinit = field_init_mont_fp;
    }
}

int pbc_mpz_set_str(mpz_t z, char *s, int base)
{
    int b, i = 0;
    mpz_set_ui(z, 0);
    if (!base) {
	b = 10;
    } else {
	if (base < 2) return 0;
	if (base > 36) return 0;
	b = base;
    }
    for (;;) {
	int j;
	char c = s[i];
	if (!c) break;
	if (isspace(c)) {
	    i++;
	    continue;
	}
	if (isdigit(c)) {
	    j = c - '0';
	} else if (c >= 'A' && c <= 'Z') {
	    j = c - 'A';
	} else if (c >= 'a' && c <= 'z') {
	    j = c - 'a';
	} else break;

	if (j >= b) break;

	mpz_mul_ui(z, z, b);
	mpz_add_ui(z, z, j);
	i++;
    }
    return i;
}
