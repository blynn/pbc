/* $GF(3^m)     = GF(3)[x]/(x^m + x^t + 2)$
   $GF(3^{2*m}) = GF(3^m)[x]/(x^2 + 1)$
   $GF(3^{3*m}) = GF(3^m)[x]/(x^3 - x -1)$
   $GF(3^{6*m}) = GF(3^{2*m})[x]/(x^3 - x -1)$

   The "gf3_*" functions are for $GF(3)$.
   The "gf3m_*" functions are for $GF(3^m)$.
   The "gf32m_*" functions are for $GF(3^{2*m})$.
   The "gf33m_*" functions are for $GF(3^{3*m})$ and $GF(3^{6*m})$.

   (gf3m field_t).data is a pointer of struct params
   (gf3m element_t).data is a pointer of unsigned long
   (gf32m element_t).data is gf32m_ptr
   (gf33m element_t).data is gf33m_ptr */

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdint.h>
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_memory.h"
#include "pbc_field.h"

typedef unsigned long gf3;

typedef struct { /* private data of $GF(3^m)$ */
    unsigned int len; /* the number of native machine integers required to represent one GF(3^m) element */
    unsigned int m; /* the irreducible polynomial is $x^m + x^t + 2$ */
    unsigned int t; /* the irreducible polynomial is $x^m + x^t + 2$ */
    element_ptr p; /* $p$ is the irreducible polynomial. */
} params;

typedef struct {
    element_t _0, _1;
} gf32m_s;

typedef gf32m_s *gf32m_ptr;

typedef struct {
    element_t _0, _1, _2;
} gf33m_s;

typedef gf33m_s *gf33m_ptr;

#define W (sizeof(unsigned long)*8) /* number of GF(3) elements in one processor integer */
#define PARAM(e) ((params *)e->field->data)
#define LEN(e) (PARAM(e)->len)
#define SIZE(e) (LEN(e) * 2 * sizeof(unsigned long))
#define DATA1(e) ((unsigned long*)e->data)
#define DATA2(e) ((unsigned long*)e->data + LEN(e))
#define GF32M(e) ((gf32m_s *)e->data)
#define GF33M(e) ((gf33m_s *)e->data)
#define BASE(e) ((field_ptr)e->field->data)
#define print(e) {printf(#e": "); element_out_str(stdout, 10, e); printf("\n");}

static size_t gf3m_out_str(FILE *stream, int base, element_t e) {
    if (base != 10 && base != 16)
        pbc_die("only support base 10 and base 16");
    size_t size = 0;
    unsigned i;
    unsigned long *d = DATA1(e);
    for (i = 0; i < LEN(e) * 2; i++) {
        if (base == 16)
            size += fprintf(stream, "0x%lx,", d[i]);
        else
            size += fprintf(stream, "%lu,", d[i]);
    }
    return size;
}

/* $a <- 0$ */
static void gf3m_zero(element_t a) {
    memset(a->data, 0, SIZE(a));
}

static void gf3m_init(element_t e) {
    e->data = pbc_malloc(SIZE(e));
    gf3m_zero(e);
}

static void gf3m_clear(element_t e) {
    pbc_free(e->data);
}

/* $e <- a$ */
static void gf3m_assign(element_t e, element_t a) {
    memcpy(e->data, a->data, SIZE(a));
}

/* $a <- a/x$. $len$ is the number of elements in $a$ */
static void shift_down(unsigned int len, unsigned long a[]) {
    unsigned long h = 0;
    const unsigned long x = 1ul << (W - 1);
    int i;
    for (i = len - 1; i >= 0; i--) {
        unsigned long l = a[i] & 1;
        a[i] >>= 1;
        if (h)
            a[i] |= x;
        h = l;
    }
}

/* $e <- e/x$ */
static void gf3m_shift_down(element_t e) {
    shift_down(LEN(e), DATA1(e));
    shift_down(LEN(e), DATA2(e));
}

/* $a <- a*x$. $len$ is the number of elements in $a$ */
static void shift_up(unsigned int len, unsigned long a[]) {
    unsigned long l = 0;
    const unsigned long x = 1ul << (W - 1), y = x - 1;
    unsigned i;
    for (i = 0; i < len; i++) {
        unsigned long h = a[i] & x;
        a[i] = ((a[i] & y) << 1) | l;
        l = h ? 1 : 0;
    }
}

/* $e <- e*x$ */
static void gf3m_shift_up(element_t e) {
    shift_up(LEN(e), DATA1(e));
    shift_up(LEN(e), DATA2(e));
}

/* return the coefficient of $x^pos$ in $e$ */
static unsigned gf3m_get(element_t e, unsigned pos) {
    unsigned long *a1 = DATA1(e), *a2 = DATA2(e);
    unsigned x = pos / W;
    unsigned long y = 1ul << (pos % W), v1 = a1[x] & y, v2 = a2[x] & y;
    return v1 ? 1 : (v2 ? 2 : 0);
}

/* set the coefficient of $x^pos$ as 1 */
static void gf3m_set(element_t e, unsigned pos, unsigned value) {
    unsigned long *a = DATA1(e);
    /* assert value == 0, 1 or 2 */
    if (value == 2)
        a = DATA2(e);
    if (value)
        a[pos / W] |= 1ul << (pos % W);
}

/* $e <- a+b$ */
static void gf3m_add(element_t e, element_t a, element_t b) {
    unsigned long *e1 = DATA1(e), *e2 = DATA2(e), *a1 = DATA1(a),
            *a2 = DATA2(a), *b1 = DATA1(b), *b2 = DATA2(b);
    unsigned i;
    for (i = 0; i < LEN(e); i++, e1++, e2++, a1++, a2++, b1++, b2++) {
        unsigned long t = (*a1 | *a2) & (*b1 | *b2), c1 = t ^ (*a1 | *b1), c2 =
                t ^ (*a2 | *b2);
        *e1 = c1;
        *e2 = c2;
    }
}

/* $e <- x-y$ */
static void gf3m_sub(element_t e, element_t a, element_t b) {
    unsigned long *e1 = DATA1(e), *e2 = DATA2(e), *a1 = DATA1(a),
            *a2 = DATA2(a), *b1 = DATA2(b), *b2 = DATA1(b);
    unsigned i;
    for (i = 0; i < LEN(e); i++, e1++, e2++, a1++, a2++, b1++, b2++) {
        unsigned long t = (*a1 | *a2) & (*b1 | *b2), c1 = t ^ (*a1 | *b1), c2 =
                t ^ (*a2 | *b2);
        *e1 = c1;
        *e2 = c2;
    }
}

/* return 0 if $a == b$ in $GF(3^m)$, 1 otherwise. */
static int gf3m_cmp(element_t a, element_t b) {
    unsigned long *pa = DATA1(a), *pb = DATA1(b);
    unsigned i;
    for (i = 0; i < LEN(a) * 2; i++, pa++, pb++)
        if (*pa != *pb)
            return 1;
    return 0;
}

/* $a <- 1$ */
static void gf3m_one(element_t a) {
    gf3m_zero(a);
    *DATA1(a) = 1;
}

static int gf3m_is0(element_t e) {
    unsigned i;
    for (i = 0; i < LEN(e) * 2; i++)
        if (DATA1(e)[i])
            return 0;
    return 1;
}

static int gf3m_is1(element_t e) {
    unsigned i;
    if (DATA1(e)[0] != 1)
        return 0;
    for (i = 1; i < LEN(e) * 2; i++)
        if (DATA1(e)[i])
            return 0;
    return 1;
}

/* set $a$ to be a random element in $GF(3^m)$ */
static void gf3m_random(element_t a) {
    /* TODO: use uniform distribution? */
    params *c = PARAM(a);
    unsigned rm = c->m % W;
    const unsigned long i1 = ~0ul;
    unsigned long i2 = (1ul << rm) - 1;
    unsigned long *a1 = DATA1(a), *a2 = DATA2(a);
    unsigned i;
    for (i = 0; i < c->len - 1; i++, a1++, a2++) { /* TODO: if $RAND_MAX < i1$ ? */
        *a1 = rand() & i1;
        *a2 = rand() & i1 & ~(*a1); /* assuring there is no bit that a1[x] & a2[x] == 1 */
    }
    unsigned long x = rm ? i2 : i1;
    *a1 = rand() & x;
    *a2 = rand() & x & ~(*a1);
}

static void swap(unsigned long *a, unsigned long *b) {
    *a ^= *b;
    *b ^= *a;
    *a ^= *b;
}

/* $y <- (-x)$ */
static void gf3m_neg(element_t y, element_t x) {
    unsigned long *a1 = DATA1(x), *a2 = DATA2(x), *c1 = DATA1(y),
            *c2 = DATA2(y);
    if (a1 == c1) {
        unsigned i;
        for (i = 0; i < LEN(y); i++, a1++, a2++)
            swap(a1, a2);
    } else {
        memcpy(c1, a2, SIZE(y) / 2);
        memcpy(c2, a1, SIZE(y) / 2);
    }
}

/* doing reduction
 * The function returns the value of $a$ modulo $the irreducible trinomial$.
 * $degree$ equals the degree of $a$.
 * $2*len$ is the number of elements in $a$ */
static void gf3m_reduct(element_t e, unsigned len, unsigned degree) {
    // the $len$ argument exists because sometimes $len != p->len$
    params *p = PARAM(e);
    unsigned old = p->len;
    p->len = len;
    element_t px;
    element_init(px, e->field);
    gf3m_set(px, degree, 1);
    gf3m_set(px, degree - p->m + p->t, 1);
    gf3m_set(px, degree - p->m, 2);
    while (degree >= p->m) {
        unsigned v = gf3m_get(e, degree);
        if (v == 1)
            gf3m_sub(e, e, px);
        else if (v == 2)
            gf3m_add(e, e, px);
        degree--;
        gf3m_shift_down(px);
    }
    element_clear(px);
    p->len = old;
}

/* doing multiplication of $n \in \{0,1,2\}$ and $a$ in $GF(3^m)$
 * The function sets $e <- n * a$. */
static void gf3m_f1(element_t e, unsigned n, element_t a) {
    /* assert $e$ is not $a$ */
    if (n == 0)
        memset(DATA1(e), 0, SIZE(e));
    else if (n == 1)
        memcpy(DATA1(e), DATA1(a), SIZE(e));
    else {
        memcpy(DATA1(e), DATA2(a), SIZE(e) / 2);
        memcpy(DATA2(e), DATA1(a), SIZE(e) / 2);
    }
}

/* $e <- e*x mod p(x)$ */
static void gf3m_f2(element_t e) {
    params *p = PARAM(e);
    gf3m_shift_up(e);
    unsigned v = gf3m_get(e, p->m);
    if (v == 1)
        gf3m_sub(e, e, p->p);
    else if (v == 2)
        gf3m_add(e, e, p->p);
}

/* doing multiplication in GF(3^m)
 * The function sets $e == a*b \in GF(3^m)$ */
static void gf3m_mult(element_t e, element_ptr a, element_t b) {
    params *p = PARAM(a);
    element_t aa, t, c;
    element_init(aa, a->field);
    element_set(aa, a);
    a = aa; // clone $a$
    element_init(t, a->field);
    element_init(c, a->field);
    unsigned i;
    for (i = 0; i < p->m; i++) {
        unsigned v = gf3m_get(b, i);
        gf3m_f1(t, v, a); /* t == b[i]*a in GF(3^m) */
        gf3m_add(c, c, t); /* c += b[i]*a in GF(3^m) */
        gf3m_f2(a); /* a == a*x in GF(3^m) */
    }
    element_set(e, c);
    element_clear(t);
    element_clear(c);
    element_clear(aa);
}

/* $e <- x^3$ */
static void gf3m_cubic(element_t e, element_t x) {
    /* TODO: faster algorithm */
    params *p = PARAM(x);
    unsigned old = p->len;
    unsigned len = (3 * p->m - 2 + W - 1) / W; /* length of $b1 */
    p->len = len;
    element_t a;
    element_init(a, x->field);
    unsigned i;
    for (i = 0; i < p->m; i++) {
        p->len = old;
        unsigned v = gf3m_get(x, i);
        p->len = len;
        gf3m_set(a, 3 * i, v);
    }
    gf3m_reduct(a, len, 3 * p->m - 3);
    p->len = old;
    memcpy(DATA1(e), DATA1(a), SIZE(e) / 2);
    memcpy(DATA2(e), DATA1(a) + len, SIZE(e) / 2);
    element_clear(a);
}

/* multiplication modulo 3 of two elements in GF(3)
 * for example, $mult(2,2) == 1$, and $mult(1,2) == 2$ */
static unsigned gf3_mult(unsigned a, unsigned b) {
    static const unsigned l[] = { 0, 1, 2, 0, 1 };
    return l[a * b];
}

static void gf3m_swap(element_t a, element_t b) {
    unsigned long *p = DATA1(a);
    a->data = b->data;
    b->data = p;
}

/* computing the inversion of an element $a$ in GF(3^m), i.e., $e <- a^{-1}$
 The algorithm is by Tim Kerins, Emanuel Popovici and William Marnane
 in the paper of "Algorithms and Architectures for use in FPGA",
 Lecture Notes in Computer Science, 2004, Volume 3203/2004, 74-83.
 Note that $U$ must have an extra bit, i.e, (_m + W - 1) // W == (_m + W) // W */
static void gf3m_invert(element_t e, element_t a) {
    struct field_s *f = a->field;
    params *p = PARAM(a);
    unsigned lenA = p->len;
    unsigned lenS = (3 * p->m + W - 1) / W;
    p->len = lenS;
    element_t S, R, t, U, V, t2;
    element_init(S, f);
    element_init(R, f);
    element_init(t, f);
    memcpy(DATA1(S), DATA1(p->p), lenA * sizeof(unsigned long)); /* S = p(x) */
    memcpy(DATA1(S) + lenS, DATA1(p->p) + lenA, lenA * sizeof(unsigned long));
    memcpy(DATA1(R), DATA1(a), lenA * sizeof(unsigned long)); /* R = _clone(a) */
    memcpy(DATA1(R) + lenS, DATA1(a) + lenA, lenA * sizeof(unsigned long));
    p->len = lenA;
    element_init(U, f);
    gf3m_one(U);
    element_init(V, f);
    element_init(t2, f);
    unsigned d = 0, i, r_m, s_m, q, x;
    for (i = 0; i < p->m * 2; i++) {
        p->len = lenS;
        r_m = gf3m_get(R, p->m), s_m = gf3m_get(S, p->m);
        if (r_m == 0) {
            gf3m_shift_up(R); /* R = xR */
            p->len = lenA;
            gf3m_f2(U); /* U = xU mod p */
            d++;
        } else {
            q = gf3_mult(r_m, s_m);
            gf3m_f1(t, q, R);
            gf3m_sub(S, S, t); /* S = S-qR */
            gf3m_shift_up(S); /* S = xS */
            p->len = lenA;
            gf3m_f1(t2, q, U);
            gf3m_sub(V, V, t2); /* V = V-qU */
            if (d == 0) {
                gf3m_swap(S, R);
                gf3m_swap(U, V);
                gf3m_f2(U); /* U = xU mod p*/
                d++;
            } else {
                x = gf3m_get(U, 0);
                if (x == 1) /* assuring x|U */
                    gf3m_add(U, U, p->p);
                else if (x == 2)
                    gf3m_sub(U, U, p->p);
                gf3m_shift_down(U); /* divide U by $x$ */
                d--;
            }
        }
    }
    p->len = lenS;
    r_m = gf3m_get(R, p->m); /* assume r_m is not zero */
    p->len = lenA;
    if (r_m == 2)
        gf3m_neg(U, U);
    memcpy(e->data, U->data, lenA * 2 * sizeof(unsigned long));
    element_clear(S);
    element_clear(R);
    element_clear(U);
    element_clear(V);
    element_clear(t);
    element_clear(t2);
}

static void gf3m_sqrt(element_t e, element_t a) {
    field_ptr f = e->field;
    mpz_t t;
    mpz_init(t); // t == (field_order  + 1) / 4
    mpz_set(t, f->order);
    mpz_add_ui(t, t, 1);
    mpz_tdiv_q_2exp(t, t, 2);
    element_pow_mpz(e, a, t);
    mpz_clear(t);
}

int gf3m_to_bytes(unsigned char *d, element_ptr e) {
    unsigned long *a = DATA1(e), *b = DATA2(e);
    unsigned long i, j;
    for (i = 0; i < LEN(e); i++, a++, b++) {
        for (j = 0; j < sizeof(unsigned long) * 8; j += 8) {
            *(d++) = (unsigned char) ((*a) >> j);
            *(d++) = (unsigned char) ((*b) >> j);
        }
    }
    return SIZE(e);
}

int gf3m_from_bytes(element_ptr e, unsigned char *d) {
    unsigned long *a = DATA1(e), *b = DATA2(e);
    unsigned i;
    int j;
    for (i = 0; i < LEN(e); i++, a++, b++, d += sizeof(unsigned long) * 2) {
        *a = 0, *b = 0;
        j = 2 * sizeof(unsigned long) - 2;
        while (j >= 0) {
            *a <<= 8, *b <<= 8;
            *a += d[j];
            *b += d[j + 1];
            j -= 2;
        }
    }
    return SIZE(e);
}

static void field_clear_gf3m(field_t f) {
    params *p = f->data;
    gf3m_clear(p->p);
    pbc_free(p->p);
    pbc_free(p);
}

/* initialize the finite field as $GF(3^m)$, whose irreducible polynomial is with the degree of $m$ */
void field_init_gf3m(field_t f, unsigned m, unsigned t) {
    params *p = pbc_malloc(sizeof(*p));
    p->len = (m + (W - 1) + 1) / W; /* extra one bit for $_p$ */
    p->m = m;
    p->t = t;
    p->p = pbc_malloc(sizeof(*(p->p)));
    p->p->field = f;
    p->p->data = pbc_malloc(2 * sizeof(unsigned long) * p->len);
    memset(p->p->data, 0, 2 * sizeof(unsigned long) * p->len);
    unsigned long *p1 = p->p->data, *p2 = p1 + p->len;
    p2[0] = 1; /* _p == x^m+x^t+2 */
    unsigned int p_t = p->t;
    p1[p_t / W] |= 1ul << (p_t % W);
    p1[m / W] |= 1ul << (m % W);

    field_init(f);
    f->field_clear = field_clear_gf3m;
    f->init = gf3m_init;
    f->clear = gf3m_clear;
    f->set = gf3m_assign;
    f->set0 = gf3m_zero;
    f->set1 = gf3m_one;
    f->is0 = gf3m_is0;
    f->is1 = gf3m_is1;
    f->add = gf3m_add;
    f->sub = gf3m_sub;
    f->mul = gf3m_mult;
    f->cubic = gf3m_cubic;
    f->invert = gf3m_invert;
    f->neg = gf3m_neg;
    f->random = gf3m_random;
    f->cmp = gf3m_cmp;
    f->sqrt = gf3m_sqrt;
    f->from_bytes = gf3m_from_bytes;
    f->to_bytes = gf3m_to_bytes;
    f->out_str = gf3m_out_str;
    f->fixed_length_in_bytes = 2 * sizeof(unsigned long) * p->len;
    f->data = p;
    f->name = "GF(3^m)";

    mpz_set_ui(f->order, 3);
    mpz_pow_ui(f->order, f->order, p->m);
}

static size_t gf32m_out_str(FILE *stream, int base, element_t e) {
    UNUSED_VAR(base);
    element_ptr e0 = GF32M(e)->_0, e1 = GF32M(e)->_1;
    size_t size = 0;
    size += element_out_str(stream, base, e0);
    size += element_out_str(stream, base, e1);
    return size;
}

static void gf32m_init(element_t e) {
    e->data = pbc_malloc(sizeof(gf32m_s));
    gf32m_ptr p = (gf32m_ptr) e->data;
    field_ptr base = BASE(e);
    element_init(p->_0, base);
    element_init(p->_1, base);
}

static void gf32m_clear(element_t e) {
    gf32m_ptr p = (gf32m_ptr) e->data;
    element_clear(p->_0);
    element_clear(p->_1);
    pbc_free(e->data);
}

static void gf32m_set0(element_t e) {
    element_ptr e0 = GF32M(e)->_0, e1 = GF32M(e)->_1;
    element_set0(e0);
    element_set0(e1);
}

static void gf32m_set1(element_t e) {
    element_ptr e0 = GF32M(e)->_0, e1 = GF32M(e)->_1;
    element_set1(e0);
    element_set0(e1);
}

static int gf32m_item_count(element_t e) {
    UNUSED_VAR(e);
    return 2;
}

static element_ptr gf32m_item(element_t a, int i) {
    element_ptr a0 = GF32M(a)->_0, a1 = GF32M(a)->_1;
    return i == 0 ? a0 : a1;
}

static void gf32m_assign(element_t e, element_t a) {
    element_ptr a0 = GF32M(a)->_0, a1 = GF32M(a)->_1, e0 = GF32M(e)->_0, e1 =
            GF32M(e)->_1;
    element_set(e0, a0);
    element_set(e1, a1);
}

static void gf32m_random(element_t e) {
    element_ptr e0 = GF32M(e)->_0, e1 = GF32M(e)->_1;
    element_random(e0);
    element_random(e1);
}

/* return 0 if $a == b$, 1 otherwise */
static int gf32m_cmp(element_t a, element_t b) {
    element_ptr a0 = GF32M(a)->_0, a1 = GF32M(a)->_1, b0 = GF32M(b)->_0, b1 =
            GF32M(b)->_1;
    return element_cmp(a0, b0) || element_cmp(a1, b1);
}

/* $c <- a+b$ */
static void gf32m_add(element_t c, element_t a, element_t b) {
    element_ptr a0 = GF32M(a)->_0, a1 = GF32M(a)->_1, b0 = GF32M(b)->_0, b1 =
            GF32M(b)->_1, c0 = GF32M(c)->_0, c1 = GF32M(c)->_1;
    element_add(c0, a0, b0);
    element_add(c1, a1, b1);
}

/* $c <- a-b$ */
static void gf32m_sub(element_t c, element_t a, element_t b) {
    element_ptr a0 = GF32M(a)->_0, a1 = GF32M(a)->_1, b0 = GF32M(b)->_0, b1 =
            GF32M(b)->_1, c0 = GF32M(c)->_0, c1 = GF32M(c)->_1;
    element_sub(c0, a0, b0);
    element_sub(c1, a1, b1);
}

/* $c <- (-a)$ */
static void gf32m_neg(element_t c, element_t a) {
    element_ptr a0 = GF32M(a)->_0, a1 = GF32M(a)->_1, c0 = GF32M(c)->_0, c1 =
            GF32M(c)->_1;
    element_neg(c0, a0);
    element_neg(c1, a1);
}

/* $e<- a*b$ */
static void gf32m_mult(element_t e, element_t a, element_t b) {
    element_ptr a0 = GF32M(a)->_0, a1 = GF32M(a)->_1, b0 = GF32M(b)->_0, b1 =
            GF32M(b)->_1, e0 = GF32M(e)->_0, e1 = GF32M(e)->_1;
    field_ptr base = BASE(a);
    element_t a0b0, a1b1, t0, t1, c1;
    element_init(a0b0, base);
    element_init(a1b1, base);
    element_init(t0, base);
    element_init(t1, base);
    element_init(c1, base);
    element_mul(a0b0, a0, b0);
    element_mul(a1b1, a1, b1);
    element_add(t0, a1, a0);
    element_add(t1, b1, b0);
    element_mul(c1, t0, t1); // c1 == (a1+a0)*(b1+b0)
    element_sub(c1, c1, a1b1);
    element_sub(c1, c1, a0b0);
    element_ptr c0 = a0b0;
    element_sub(c0, c0, a1b1); // c0 == a0*b0 - a1*b1
    element_set(e0, c0);
    element_set(e1, c1);
    element_clear(a0b0);
    element_clear(a1b1);
    element_clear(t0);
    element_clear(t1);
    element_clear(c1);
}

/* $e <- a^3$ */
static void gf32m_cubic(element_t e, element_t a) {
    element_ptr a0 = GF32M(a)->_0, a1 = GF32M(a)->_1, e0 = GF32M(e)->_0, e1 =
            GF32M(e)->_1;
    field_ptr base = BASE(a);
    element_t c0, c1;
    element_init(c0, base);
    element_init(c1, base);
    element_cubic(c0, a0);
    element_cubic(c1, a1);
    element_neg(c1, c1); // c1 == -(a1^3)
    element_set(e0, c0);
    element_set(e1, c1);
    element_clear(c0);
    element_clear(c1);
}

void field_clear_gf32m(field_t f) {
    UNUSED_VAR(f);
}

/* initialize the finite field as $base_field[x]/(x^2 + 1)$, whose base field is $b$ */
void field_init_gf32m(field_t f, field_t b) {
    field_init(f);
    f->data = b;
    f->field_clear = field_clear_gf32m;
    f->init = gf32m_init;
    f->clear = gf32m_clear;
    f->set = gf32m_assign;
    f->set0 = gf32m_set0;
    f->set1 = gf32m_set1;
    f->random = gf32m_random;
    f->cmp = gf32m_cmp;
    f->add = gf32m_add;
    f->sub = gf32m_sub;
    f->neg = gf32m_neg;
    f->mul = gf32m_mult;
    f->cubic = gf32m_cubic;
    f->item_count = gf32m_item_count;
    f->item = gf32m_item;
    f->out_str = gf32m_out_str;
    mpz_pow_ui(f->order, b->order, 2);
    f->name = "GF(3^{2*m})";
}

static size_t gf33m_out_str(FILE *stream, int base, element_t e) {
    UNUSED_VAR(base);
    element_ptr e0 = GF33M(e)->_0, e1 = GF33M(e)->_1, e2 = GF33M(e)->_2;
    size_t size = 0;
    size += element_out_str(stream, base, e0);
    size += element_out_str(stream, base, e1);
    size += element_out_str(stream, base, e2);
    return size;
}

static void gf33m_init(element_t e) {
    e->data = pbc_malloc(sizeof(gf33m_s));
    gf33m_ptr p = (gf33m_ptr) e->data;
    field_ptr base = BASE(e);
    element_init(p->_0, base);
    element_init(p->_1, base);
    element_init(p->_2, base);
}

static void gf33m_clear(element_t e) {
    gf33m_ptr p = (gf33m_ptr) e->data;
    element_clear(p->_0);
    element_clear(p->_1);
    element_clear(p->_2);
    pbc_free(e->data);
}

static void gf33m_set0(element_t e) {
    element_ptr e0 = GF33M(e)->_0, e1 = GF33M(e)->_1, e2 = GF33M(e)->_2;
    element_set0(e0);
    element_set0(e1);
    element_set0(e2);
}

static void gf33m_set1(element_t e) {
    element_ptr e0 = GF33M(e)->_0, e1 = GF33M(e)->_1, e2 = GF33M(e)->_2;
    element_set1(e0);
    element_set0(e1);
    element_set0(e2);
}

static int gf33m_item_count(element_t e) {
    UNUSED_VAR(e);
    return 3;
}

static element_ptr gf33m_item(element_t a, int i) {
    element_ptr a0 = GF33M(a)->_0, a1 = GF33M(a)->_1, a2 = GF33M(a)->_2;
    return i == 0 ? a0 : (i == 1 ? a1 : a2);
}

static void gf33m_assign(element_t e, element_t a) {
    element_ptr a0 = GF33M(a)->_0, a1 = GF33M(a)->_1, a2 = GF33M(a)->_2, e0 =
            GF33M(e)->_0, e1 = GF33M(e)->_1, e2 = GF33M(e)->_2;
    element_set(e0, a0);
    element_set(e1, a1);
    element_set(e2, a2);
}

static void gf33m_random(element_t e) {
    element_ptr e0 = GF33M(e)->_0, e1 = GF33M(e)->_1, e2 = GF33M(e)->_2;
    element_random(e0);
    element_random(e1);
    element_random(e2);
}

/* return 0 if $a == b$, 1 otherwise */
static int gf33m_cmp(element_t a, element_t b) {
    element_ptr a0 = GF33M(a)->_0, a1 = GF33M(a)->_1, a2 = GF33M(a)->_2, b0 =
            GF33M(b)->_0, b1 = GF33M(b)->_1, b2 = GF33M(b)->_2;
    return element_cmp(a0, b0) || element_cmp(a1, b1) || element_cmp(a2, b2);
}

/* $c <- a+b$ */
static void gf33m_add(element_t c, element_t a, element_t b) {
    element_ptr a0 = GF33M(a)->_0, a1 = GF33M(a)->_1, a2 = GF33M(a)->_2, b0 =
            GF33M(b)->_0, b1 = GF33M(b)->_1, b2 = GF33M(b)->_2, c0 =
            GF33M(c)->_0, c1 = GF33M(c)->_1, c2 = GF33M(c)->_2;
    element_add(c0, a0, b0);
    element_add(c1, a1, b1);
    element_add(c2, a2, b2);
}

/* $c <- a-b$ */
static void gf33m_sub(element_t c, element_t a, element_t b) {
    element_ptr a0 = GF33M(a)->_0, a1 = GF33M(a)->_1, a2 = GF33M(a)->_2, b0 =
            GF33M(b)->_0, b1 = GF33M(b)->_1, b2 = GF33M(b)->_2, c0 =
            GF33M(c)->_0, c1 = GF33M(c)->_1, c2 = GF33M(c)->_2;
    element_sub(c0, a0, b0);
    element_sub(c1, a1, b1);
    element_sub(c2, a2, b2);
}

/* $c <- a*b$ */
static void gf33m_mult(element_t e, element_t a, element_t b) {
    element_ptr a0 = GF33M(a)->_0, a1 = GF33M(a)->_1, a2 = GF33M(a)->_2, b0 =
            GF33M(b)->_0, b1 = GF33M(b)->_1, b2 = GF33M(b)->_2, e0 =
            GF33M(e)->_0, e1 = GF33M(e)->_1, e2 = GF33M(e)->_2;
    field_ptr base = BASE(e);
    element_t t0, t1, c1, a0b0, a1b1, a2b2;
    element_init(t0, base);
    element_init(t1, base);
    element_init(c1, base);
    element_init(a0b0, base);
    element_init(a1b1, base);
    element_init(a2b2, base);
    element_mul(a0b0, a0, b0);
    element_mul(a1b1, a1, b1);
    element_mul(a2b2, a2, b2);
    element_ptr d0 = a0b0;
    element_add(t0, a1, a0);
    element_add(t1, b1, b0);
    element_t d1;
    element_init(d1, base);
    element_mul(d1, t0, t1);
    element_sub(d1, d1, a1b1);
    element_sub(d1, d1, a0b0);
    element_add(t0, a2, a0);
    element_add(t1, b2, b0);
    element_t d2;
    element_init(d2, base);
    element_mul(d2, t0, t1);
    element_add(d2, d2, a1b1);
    element_sub(d2, d2, a2b2);
    element_sub(d2, d2, a0b0);
    element_add(t0, a2, a1);
    element_add(t1, b2, b1);
    element_t d3;
    element_init(d3, base);
    element_mul(d3, t0, t1);
    element_sub(d3, d3, a2b2);
    element_sub(d3, d3, a1b1);
    element_ptr d4 = a2b2;
    element_add(t0, d0, d3);
    element_ptr c0 = t0;
    element_add(c1, d1, d3);
    element_add(c1, c1, d4);
    element_add(t1, d2, d4);
    element_ptr c2 = t1;
    element_set(e0, c0);
    element_set(e1, c1);
    element_set(e2, c2);
    element_clear(t0);
    element_clear(t1);
    element_clear(c1);
    element_clear(a0b0);
    element_clear(a1b1);
    element_clear(a2b2);
    element_clear(d1);
    element_clear(d2);
    element_clear(d3);
}

/* $e <- a^3$ */
static void gf33m_cubic(element_t e, element_t a) {
    field_ptr base = BASE(a);
    element_ptr a0 = GF33M(a)->_0, a1 = GF33M(a)->_1, a2 = GF33M(a)->_2, e0 =
            GF33M(e)->_0, e1 = GF33M(e)->_1, e2 = GF33M(e)->_2;
    element_t a03, a13, a23;
    element_init(a03, base);
    element_init(a13, base);
    element_init(a23, base);
    element_cubic(a03, a0);
    element_cubic(a13, a1);
    element_cubic(a23, a2);
    element_add(a03, a03, a13);
    element_add(a03, a03, a23);
    element_ptr c0 = a03;
    element_sub(a13, a13, a23);
    element_ptr c1 = a13;
    element_ptr c2 = a23;
    element_set(e0, c0);
    element_set(e1, c1);
    element_set(e2, c2);
    element_clear(a03);
    element_clear(a13);
    element_clear(a23);
}

/* $e <- a^{-1}$ */
static void gf33m_invert(element_t e, element_t a) {
    element_ptr a0 = GF33M(a)->_0, a1 = GF33M(a)->_1, a2 = GF33M(a)->_2, e0 =
            GF33M(e)->_0, e1 = GF33M(e)->_1, e2 = GF33M(e)->_2;
    field_ptr base = BASE(e);
    element_t a02, a12, a22;
    element_init(a02, base);
    element_init(a12, base);
    element_init(a22, base);
    element_mul(a02, a0, a0);
    element_mul(a12, a1, a1);
    element_mul(a22, a2, a2);
    element_t v0;
    element_init(v0, base);
    element_sub(v0, a0, a2); // v0 == a0-a2
    element_t delta;
    element_init(delta, base);
    element_mul(delta, v0, a02); // delta = (a0-a2)*(a0^2), free
    element_sub(v0, a1, a0); // v0 == a1-a0
    element_t c0;
    element_init(c0, base);
    element_mul(c0, v0, a12); // c0 == (a1-a0)*(a1^2)
    element_add(delta, delta, c0); // delta = (a0-a2)*(a0^2) + (a1-a0)*(a1^2)
    element_sub(v0, a2, v0); // v0 == a2-(a1-a0) = a0-a1+a2
    element_t c1;
    element_init(c1, base);
    element_mul(c1, v0, a22); // c1 == (a0-a1+a2)*(a2^2)
    element_add(delta, delta, c1); // delta = (a0-a2)*(a0^2) + (a1-a0)*(a1^2) + (a0-a1+a2)*(a2^2)
    element_invert(delta, delta); // delta = [(a0-a2)*(a0^2) + (a1-a0)*(a1^2) + (a0-a1+a2)*(a2^2)] ^ {-1}
    element_add(v0, a02, a22); // v0 == a0^2+a2^2
    element_t c2;
    element_init(c2, base);
    element_mul(c2, a0, a2); // c2 == a0*a2
    element_sub(c0, v0, c2); // c0 == a0^2+a2^2-a0*a2
    element_add(v0, a1, a2); // v0 == a1+a2
    element_t c3;
    element_init(c3, base);
    element_mul(c3, a1, v0); // c3 == a1*(a1+a2)
    element_sub(c0, c0, c3); // c0 == a0^2+a2^2-a0*a2-a1*(a1+a2)
    element_mul(c0, c0, delta); // c0 *= delta
    element_mul(c1, a0, a1); // c1 == a0*a1
    element_sub(c1, a22, c1); // c1 == a2^2-a0*a1
    element_mul(c1, c1, delta); // c1 *= delta
    element_sub(c2, a12, c2); // c2 == a1^2-a0*a2
    element_sub(c2, c2, a22); // c2 == a1^2-a0*a2-a2^2
    element_mul(c2, c2, delta); // c2 *= delta
    element_set(e0, c0);
    element_set(e1, c1);
    element_set(e2, c2);
    element_clear(a02);
    element_clear(a12);
    element_clear(a22);
    element_clear(v0);
    element_clear(delta);
    element_clear(c0);
    element_clear(c1);
    element_clear(c2);
    element_clear(c3);
}

void field_clear_gf33m(field_t f) {
    UNUSED_VAR(f);
}

/* initialize the finite field as $base_field[x]/(x^3 - x - 1)$, whose base field is $b$ */
void field_init_gf33m(field_t f, field_t b) {
    field_init(f);
    f->data = b;
    f->field_clear = field_clear_gf33m;
    f->init = gf33m_init;
    f->clear = gf33m_clear;
    f->set = gf33m_assign;
    f->set0 = gf33m_set0;
    f->set1 = gf33m_set1;
    f->random = gf33m_random;
    f->cmp = gf33m_cmp;
    f->add = gf33m_add;
    f->sub = gf33m_sub;
    f->mul = gf33m_mult;
    f->cubic = gf33m_cubic;
    f->invert = gf33m_invert;
    f->item_count = gf33m_item_count;
    f->item = gf33m_item;
    f->out_str = gf33m_out_str;
    mpz_pow_ui(f->order, b->order, 3);
    f->name = "GF(3^{3*m})";
}

