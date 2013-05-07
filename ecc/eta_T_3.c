#include <stdarg.h>
#include <stdio.h>
#include <stdint.h>
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_field.h"
#include "pbc_fp.h"
#include "pbc_memory.h"
#include "pbc_param.h"
#include "pbc_pairing.h"
#include "pbc_ternary_extension_field.h"
#include "param.h"

typedef struct { /* private data of $GF(3^m)$ */
    unsigned int len; /* the number of native machine integers required to represent one GF(3^m) element */
    int m; /* the irreducible polynomial is $x^m + x^t + 2$ */
    int t; /* the irreducible polynomial is $x^m + x^t + 2$ */
    element_ptr p; /* $p$ is the irreducible polynomial. */
    mpz_t n; /* group order of $G_1$, $G_2$, $G_T$ */
    mpz_t n2; /* order(elliptic curve points) / order(G_1) */
} params;

struct pairing_data {
    field_t gf3m, gf32m, gf36m;
    mpz_t n2; // cofactor
};
typedef struct pairing_data *pairing_data_ptr;

#define PARAM(e) ((params *)e->field->data)
#define ITEM(e,x,y) (element_item(element_item(e,x),y))
#define print(e) {printf(#e": "); element_out_str(stdout, 10, e); printf("\n");}

struct point_s { // points on the elliptic curve $y^2=x^3-x+1$
    int isinf;
    element_t x, y;
};
typedef struct point_s *point_ptr;
typedef struct point_s point_t[1];

#define FIELD(e) ((field_ptr) e->field)
#define BASE(e) ((field_ptr) FIELD(e)->data)
#define DATA(e) ((point_ptr) e->data)

static void point_set(element_t e, element_t a) {
    point_ptr r = DATA(e), p = DATA(a);
    r->isinf = p->isinf;
    if (!p->isinf) {
        element_set(r->x, p->x);
        element_set(r->y, p->y);
    }
}

static void point_init(element_t e) {
    field_ptr f = BASE(e);
    e->data = pbc_malloc(sizeof(struct point_s));
    point_ptr p = DATA(e);
    element_init(p->x, f);
    element_init(p->y, f);
    p->isinf = 1;
}

static void point_clear(element_t e) {
    point_ptr p = DATA(e);
    element_clear(p->x);
    element_clear(p->y);
    pbc_free(p);
}

/* return 1 if $a!=b$, 0 otherwise. */
static int point_cmp(element_t a, element_t b) {
    point_ptr pa = DATA(a), pb = DATA(b);
    if (pa->isinf == pb->isinf) {
        if (pa->isinf)
            return 0;
        else
            return element_cmp(pa->x, pb->x) || element_cmp(pa->y, pb->y);
    } else
        return 1;
}

static void point_set0(element_ptr e) {
    DATA(e)->isinf = 1;
}

static int point_is0(element_ptr e) {
    return DATA(e)->isinf;
}

static void point_random(element_t a) {
    point_ptr p = DATA(a);
    element_ptr x = p->x, y = p->y;
    field_ptr f = x->field;
    p->isinf = 0;
    element_t t, t2, e1;
    element_init(t, f);
    element_init(e1, f);
    element_set1(e1);
    element_init(t2, f);
    do {
        element_random(x);
        if (element_is0(x))
            continue;
        element_cubic(t, x); // t == x^3
        element_sub(t, t, x); // t == x^3 - x
        element_add(t, t, e1); // t == x^3 - x + 1
        element_sqrt(y, t);  // y == sqrt(x^3 - x + 1)
        element_mul(t2, y, y); // t2 == x^3 - x + 1
    } while (element_cmp(t2, t)); // t2 != t

    // make sure order of $a$ is order of $G_1$
    pairing_ptr pairing = FIELD(a)->pairing;
    pairing_data_ptr dp = pairing->data;
    element_pow_mpz(a, a, dp->n2);

    element_clear(t);
    element_clear(t2);
    element_clear(e1);
}

static void point_add(element_t c, element_t a, element_t b) {
    point_ptr p1 = DATA(a), p2 = DATA(b), p3 = DATA(c);
    int inf1 = p1->isinf, inf2 = p2->isinf;
    element_ptr x1 = p1->x, y1 = p1->y, x2 = p2->x, y2 = p2->y;
    field_ptr f = FIELD(x1);
    if (inf1) {
        point_set(c, b);
        return;
    }
    if (inf2) {
        point_set(c, a);
        return;
    }
    element_t v0, v1, v2, v3, v4, ny2;
    element_init(v0, f);
    element_init(v1, f);
    element_init(v2, f);
    element_init(v3, f);
    element_init(v4, f);
    element_init(ny2, f);
    if (!element_cmp(x1, x2)) { // x1 == x2
        element_neg(ny2, y2); // ny2 == -y2
        if (!element_cmp(y1, ny2)) {
            p3->isinf = 1;
            goto end;
        }
        if (!element_cmp(y1, y2)) { // y1 == y2
            element_invert(v0, y1); // v0 == y1^{-1}
            element_mul(v1, v0, v0); // v1 == [y1^{-1}]^2
            element_add(p3->x, v1, x1); // v1 == [y1^{-1}]^2 + x1
            element_cubic(v2, v0); // v2 == [y1^{-1}]^3
            element_add(v2, v2, y1); // v2 == [y1^{-1}]^3 + y1
            element_neg(p3->y, v2); // p3 == -([y1^{-1}]^3 + y1)
            p3->isinf = 0;
            goto end;
        }
    }
    // $P1 \ne \pm P2$
    element_sub(v0, x2, x1); // v0 == x2-x1
    element_invert(v1, v0); // v1 == (x2-x1)^{-1}
    element_sub(v0, y2, y1); // v0 == y2-y1
    element_mul(v2, v0, v1); // v2 == (y2-y1)/(x2-x1)
    element_mul(v3, v2, v2); // v3 == [(y2-y1)/(x2-x1)]^2
    element_cubic(v4, v2); // v4 == [(y2-y1)/(x2-x1)]^3
    element_add(v0, x1, x2); // v0 == x1+x2
    element_sub(v3, v3, v0); // v3 == [(y2-y1)/(x2-x1)]^2 - (x1+x2)
    element_add(v0, y1, y2); // v0 == y1+y2
    element_sub(v4, v0, v4); // v4 == (y1+y2) - [(y2-y1)/(x2-x1)]^3
    p3->isinf = 0;
    element_set(p3->x, v3);
    element_set(p3->y, v4);
    end: element_clear(v0);
    element_clear(v1);
    element_clear(v2);
    element_clear(v3);
    element_clear(v4);
    element_clear(ny2);
}

static void point_invert(element_ptr e, element_ptr a) {
    point_ptr r = DATA(e), p = DATA(a);
    r->isinf = p->isinf;
    if (!p->isinf) {
        element_set(r->x, p->x);
        element_neg(r->y, p->y);
    }
}

static size_t point_out_str(FILE *stream, int base, element_ptr a) {
    point_ptr p = DATA(a);
    size_t size = 0;
    if (p->isinf)
        return fprintf(stream, "O");
    else {
        size += element_out_str(stream, base, p->x);
        size += element_out_str(stream, base, p->y);
        return size;
    }
}

static void point_field_clear(field_ptr f) {
    UNUSED_VAR(f);
}

void field_init_eta_T_3(field_t f, field_t base) {
    field_init(f);
    f->data = (void *) base;
    f->init = point_init;
    f->clear = point_clear;
    f->random = point_random;
    f->set = point_set;
    f->cmp = point_cmp;
    f->invert = f->neg = point_invert;
    f->mul = f->add = point_add;
    f->set1 = f->set0 = point_set0;
    f->is1 = f->is0 = point_is0;
    f->mul_mpz = f->pow_mpz;
    f->out_str = point_out_str;
    f->field_clear = point_field_clear;
    f->name = "eta_T_3 point group";
}

/* computing of $(-t^2 +u*s -t*p -p^2)^3$
 * The algorithm is by J.Beuchat et.al, in the paper of "Algorithms and Arithmetic Operators for Computing
 * the $eta_T$ Pairing in Characteristic Three", algorithm 4 in the appendix */
static void algorithm4a(element_t S, element_t t, element_t u) {
    field_ptr f = FIELD(t);
    element_t e1, c0, c1, m0, v0, v2;
    element_init(e1, f);
    element_init(c0, f);
    element_init(c1, f);
    element_init(m0, f);
    element_init(v0, f);
    element_init(v2, f);
    element_set1(e1);
    element_cubic(c0, t); // c0 == t^3
    element_cubic(c1, u);
    element_neg(c1, c1); // c1 == -u^3
    element_mul(m0, c0, c0); // m0 == c0^2
    element_neg(v0, m0); // v0 == -c0^2
    element_sub(v0, v0, c0); // v0 == -c0^2 -c0
    element_sub(v0, v0, e1); // v0 == -c0^2 -c0 -1
    element_set1(v2);
    element_sub(v2, v2, c0); // v2 == 1 -c0
    // v1 == c1
    // S == [[v0, v1], [v2, f3m.zero()], [f3m.two(), f3m.zero()]]
    element_set(ITEM(S,0,0), v0);
    element_set(ITEM(S,0,1), c1);
    element_set(ITEM(S,1,0), v2);
    element_set0(ITEM(S,1,1));
    element_neg(ITEM(S,2,0), e1);
    element_set0(ITEM(S,2,1));
    element_clear(e1);
    element_clear(c0);
    element_clear(c1);
    element_clear(m0);
    element_clear(v0);
    element_clear(v2);
}

static void algorithm5(element_t c, element_ptr xp, element_ptr yp,
        element_ptr xq, element_ptr yq) {
    params *p = PARAM(xp);
    unsigned int re = p->m % 12;
    field_ptr f = FIELD(xp) /*GF(3^m)*/, f6 = FIELD(c) /*GF(3^{6*m})*/;
    element_t e1, xpp, ypp, xqq, yqq, t, nt, nt2, v1, v2, a1, a2, R, u, nu, S, S2;
    element_init(e1, f);
    element_init(xpp, f);
    element_init(ypp, f);
    element_init(xqq, f);
    element_init(yqq, f);
    element_init(t, f);
    element_init(nt, f);
    element_init(nt2, f);
    element_init(v1, f);
    element_init(v2, f);
    element_init(a1, f6);
    element_init(a2, f6);
    element_init(R, f6);
    element_init(u, f);
    element_init(nu, f);
    element_init(S, f6);
    element_init(S2, f6);
    element_set1(e1);
    element_set(xpp, xp);
    xp = xpp; // clone
    element_add(xp, xp, e1); // xp == xp + b
    element_set(ypp, yp);
    yp = ypp; // clone
    if (re == 1 || re == 11)
        element_neg(yp, yp); // yp == -\mu*b*yp, \mu == 1 when re==1, or 11
    element_set(xqq, xq);
    xq = xqq; // clone
    element_cubic(xq, xq); // xq == xq^3
    element_set(yqq, yq);
    yq = yqq; // clone
    element_cubic(yq, yq); // yq == yq^3
    element_add(t, xp, xq); // t == xp+xq
    element_neg(nt, t); // nt == -t
    element_mul(nt2, t, nt); // nt2 == -t^2
    element_mul(v2, yp, yq); // v2 == yp*yq
    element_mul(v1, yp, t); // v1 == yp*t
    if (re == 7 || re == 11) { // \lambda == 1
        element_t nyp, nyq;
        element_init(nyp, f);
        element_init(nyq, f);
        element_neg(nyp, yp); // nyp == -yp
        element_neg(nyq, yq); // nyq == -yq
        element_set(ITEM(a1,0,0), v1);
        element_set(ITEM(a1,0,1), nyq);
        element_set(ITEM(a1,1,0), nyp);
        element_clear(nyp);
        element_clear(nyq);
    } else { // \lambda == -1
        element_neg(v1, v1); // v1 == -yp*t
        element_set(ITEM(a1,0,0), v1);
        element_set(ITEM(a1,0,1), yq);
        element_set(ITEM(a1,1,0), yp);
    }
    // a2 == -t^2 +yp*yq*s -t*p -p^2
    element_set(ITEM(a2,0,0), nt2);
    element_set(ITEM(a2,0,1), v2);
    element_set(ITEM(a2,1,0), nt);
    element_neg(ITEM(a2,2,0), e1);
    element_mul(R, a1, a2);
    int i;
    for (i = 0; i < (p->m - 1) / 4; i++) {
        element_cubic(R, R);
        element_cubic(R, R); // R <= R^9
        element_cubic(xq, xq);
        element_cubic(xq, xq);
        element_sub(xq, xq, e1); // xq <= xq^9-b
        element_cubic(yq, yq);
        element_cubic(yq, yq); // yq <= yq^9
        element_add(t, xp, xq); // t == xp+xq
        element_mul(u, yp, yq); // u == yp*yq
        element_neg(nu, u); // nu == -yp*yq
        algorithm4a(S, t, nu); // S == (-t^2 -u*s -t*p -p^2)^3
        element_cubic(xq, xq);
        element_cubic(xq, xq);
        element_sub(xq, xq, e1); // xq <= xq^9-b
        element_cubic(yq, yq);
        element_cubic(yq, yq); // yq <= yq^9
        element_add(t, xp, xq); // t == xp+xq
        element_mul(u, yp, yq); // u == yp*yq
        element_neg(nt, t); // nt == -t
        element_mul(nt2, t, nt); // nt2 == -t^2
        // S2 = [[nt2, u], [nt, f3m.zero()], [f3m.two(), f3m.zero()]]
        // S2 == -t^2 +u*s -t*p -p^2
        element_set(ITEM(S2,0,0), nt2);
        element_set(ITEM(S2,0,1), u);
        element_set(ITEM(S2,1,0), nt);
        element_set0(ITEM(S2,1,1));
        element_neg(ITEM(S2,2,0), e1);
        element_set0(ITEM(S2,2,1));
        element_mul(S, S, S2);
        element_mul(R, R, S);
    }
    element_set(c, R);
    element_clear(e1);
    element_clear(xpp);
    element_clear(ypp);
    element_clear(xqq);
    element_clear(yqq);
    element_clear(t);
    element_clear(nt);
    element_clear(nt2);
    element_clear(v1);
    element_clear(v2);
    element_clear(a1);
    element_clear(a2);
    element_clear(R);
    element_clear(u);
    element_clear(nu);
    element_clear(S);
    element_clear(S2);
}

/* this is the algorithm 4 in the paper of J.Beuchat et.al, "Algorithms and Arithmetic Operators for Computing
 * the $eta_T$ Pairing in Characteristic Three" */
static void algorithm4(element_t c, element_ptr xp, element_ptr yp,
        element_ptr xq, element_ptr yq) {
    params *p = PARAM(xp);
    unsigned int re = p->m % 12;
    field_ptr f = FIELD(xp) /*GF(3^m)*/, f6 = FIELD(c) /*GF(3^{6*m})*/;
    element_t e1, xpp, ypp, xqq, yqq, t, nt, nt2, v1, v2, a1, a2, R, u, S;
    element_init(e1, f);
    element_init(xpp, f);
    element_init(ypp, f);
    element_init(xqq, f);
    element_init(yqq, f);
    element_init(t, f);
    element_init(nt, f);
    element_init(nt2, f);
    element_init(v1, f);
    element_init(v2, f);
    element_init(a1, f6);
    element_init(a2, f6);
    element_init(R, f6);
    element_init(u, f);
    element_init(S, f6);
    element_set1(e1);
    element_set(xpp, xp);
    xp = xpp; // clone
    element_add(xp, xp, e1); // xp == xp + b
    element_set(ypp, yp);
    yp = ypp; // clone
    if (re == 1 || re == 11)
        element_neg(yp, yp); // yp == -\mu*b*yp, \mu == 1 when re==1, or 11
    element_set(xqq, xq);
    xq = xqq; // clone
    element_cubic(xq, xq); // xq == xq^3
    element_set(yqq, yq);
    yq = yqq; // clone
    element_cubic(yq, yq); // yq == yq^3
    element_add(t, xp, xq); // t == xp+xq
    element_neg(nt, t); // nt == -t
    element_mul(nt2, t, nt); // nt2 == -t^2
    element_mul(v2, yp, yq); // v2 == yp*yq
    element_mul(v1, yp, t); // v1 == yp*t
    if (re == 7 || re == 11) { // \lambda == 1
        element_t nyp, nyq;
        element_init(nyp, f);
        element_init(nyq, f);
        element_neg(nyp, yp); // nyp == -yp
        element_neg(nyq, yq); // nyq == -yq
        element_set(ITEM(a1,0,0), v1);
        element_set(ITEM(a1,0,1), nyq);
        element_set(ITEM(a1,1,0), nyp);
        element_clear(nyp);
        element_clear(nyq);
    } else { // \lambda == -1
        element_neg(v1, v1); // v1 == -yp*t
        element_set(ITEM(a1,0,0), v1);
        element_set(ITEM(a1,0,1), yq);
        element_set(ITEM(a1,1,0), yp);
    }
    // a2 == -t^2 +yp*yq*s -t*p -p^2
    element_set(ITEM(a2,0,0), nt2);
    element_set(ITEM(a2,0,1), v2);
    element_set(ITEM(a2,1,0), nt);
    element_neg(ITEM(a2,2,0), e1);
    element_mul(R, a1, a2);
    int i;
    for (i = 0; i < (p->m - 1) / 2; i++) {
        element_cubic(R, R);
        element_cubic(xq, xq);
        element_cubic(xq, xq);
        element_sub(xq, xq, e1); // xq <= xq^9-b
        element_cubic(yq, yq);
        element_cubic(yq, yq);
        element_neg(yq, yq); // yq <= -yq^9
        element_add(t, xp, xq); // t == xp+xq
        element_neg(nt, t); // nt == -t
        element_mul(nt2, t, nt); // nt2 == -t^2
        element_mul(u, yp, yq); // u == yp*yq
        element_set0(S);
        element_set(ITEM(S,0,0), nt2);
        element_set(ITEM(S,0,1), u);
        element_set(ITEM(S,1,0), nt);
        element_neg(ITEM(S,2,0), e1);
        element_mul(R, R, S);
    }
    element_set(c, R);
    element_clear(e1);
    element_clear(xpp);
    element_clear(ypp);
    element_clear(xqq);
    element_clear(yqq);
    element_clear(t);
    element_clear(nt);
    element_clear(nt2);
    element_clear(v1);
    element_clear(v2);
    element_clear(a1);
    element_clear(a2);
    element_clear(R);
    element_clear(u);
    element_clear(S);
}

/* computation of $c <- U ^ {3^{3m} - 1}$
 * This is the algorithm 6 in the paper above. */
static void algorithm6(element_t c, element_t u) {
    element_ptr u0 = ITEM(u,0,0), u1 = ITEM(u,0,1), u2 = ITEM(u,1,0), u3 =
            ITEM(u,1,1), u4 = ITEM(u,2,0), u5 = ITEM(u,2,1);
    field_ptr f = FIELD(u0); /*GF(3^m)*/
    field_t f3; /*GF(3^{3*m})*/
    field_init_gf33m(f3, f);
    element_t v0, v1, m0, m1, m2, a0, a1, i;
    element_init(v0, f3);
    element_init(v1, f3);
    element_init(m0, f3);
    element_init(m1, f3);
    element_init(m2, f3);
    element_init(a0, f3);
    element_init(a1, f3);
    element_init(i, f3);
    element_set(element_item(v0, 0), u0);
    element_set(element_item(v0, 1), u2);
    element_set(element_item(v0, 2), u4);
    element_set(element_item(v1, 0), u1);
    element_set(element_item(v1, 1), u3);
    element_set(element_item(v1, 2), u5);
    element_mul(m0, v0, v0);
    element_mul(m1, v1, v1);
    element_mul(m2, v0, v1);
    element_sub(a0, m0, m1);
    element_add(a1, m0, m1);
    element_invert(i, a1);
    element_mul(v0, a0, i);
    element_mul(v1, m2, i);
    element_set(ITEM(c,0,0), element_item(v0, 0));
    element_set(ITEM(c,1,0), element_item(v0, 1));
    element_set(ITEM(c,2,0), element_item(v0, 2));
    element_set(ITEM(c,0,1), element_item(v1, 0));
    element_set(ITEM(c,1,1), element_item(v1, 1));
    element_set(ITEM(c,2,1), element_item(v1, 2));
    element_clear(v0);
    element_clear(v1);
    element_clear(m0);
    element_clear(m1);
    element_clear(m2);
    element_clear(a0);
    element_clear(a1);
    element_clear(i);
    field_clear(f3);
}

/* computation of $c <- U ^ {3^m+1}$, $U \in T_2(F_{3^3M})$
 * This is the algorithm 7 in the paper above. */
static void algorithm7(element_t c, element_t u) {
    element_ptr u0 = ITEM(u,0,0), u1 = ITEM(u,0,1), u2 = ITEM(u,1,0), u3 =
            ITEM(u,1,1), u4 = ITEM(u,2,0), u5 = ITEM(u,2,1);
    field_ptr f = FIELD(u0); /*GF(3^m)*/
    params *p = PARAM(u0);
    element_t a0, a1, a2, a3, a4, a5, a6, m0, m1, m2, m3, m4, m5, m6, m7, m8,
            v0, v1, v2, v3, v4, v5, e1;
    element_init(a0, f);
    element_init(a1, f);
    element_init(a2, f);
    element_init(a3, f);
    element_init(a4, f);
    element_init(a5, f);
    element_init(a6, f);
    element_init(m0, f);
    element_init(m1, f);
    element_init(m2, f);
    element_init(m3, f);
    element_init(m4, f);
    element_init(m5, f);
    element_init(m6, f);
    element_init(m7, f);
    element_init(m8, f);
    element_init(v0, f);
    element_init(v1, f);
    element_init(v2, f);
    element_init(v3, f);
    element_init(v4, f);
    element_init(v5, f);
    element_init(e1, f);
    element_set1(e1);
    element_add(a0, u0, u1);
    element_add(a1, u2, u3);
    element_sub(a2, u4, u5);
    element_mul(m0, u0, u4);
    element_mul(m1, u1, u5);
    element_mul(m2, u2, u4);
    element_mul(m3, u3, u5);
    element_mul(m4, a0, a2);
    element_mul(m5, u1, u2);
    element_mul(m6, u0, u3);
    element_mul(m7, a0, a1);
    element_mul(m8, a1, a2);
    element_add(a3, m5, m6);
    element_sub(a3, a3, m7);
    element_neg(a4, m2);
    element_sub(a4, a4, m3);
    element_sub(a5, m3, m2);
    element_sub(a6, m1, m0);
    element_add(a6, a6, m4);
    if (p->m % 6 == 1) {
        element_add(v0, m0, m1);
        element_add(v0, v0, a4);
        element_add(v0, e1, v0);
        element_sub(v1, m5, m6);
        element_add(v1, v1, a6);
        element_sub(v2, a4, a3);
        element_add(v3, m8, a5);
        element_sub(v3, v3, a6);
        element_add(v4, a3, a4);
        element_neg(v4, v4);
        element_add(v5, m8, a5);
    } else { // p->m % 6 == 5
        element_add(v0, m0, m1);
        element_sub(v0, v0, a4);
        element_add(v0, e1, v0);
        element_sub(v1, m6, m5);
        element_add(v1, v1, a6);
        element_set(v2, a3);
        element_add(v3, m8, a5);
        element_add(v3, v3, a6);
        element_add(v4, a3, a4);
        element_neg(v4, v4);
        element_add(v5, m8, a5);
        element_neg(v5, v5);
    }
    element_set(ITEM(c,0,0), v0);
    element_set(ITEM(c,0,1), v1);
    element_set(ITEM(c,1,0), v2);
    element_set(ITEM(c,1,1), v3);
    element_set(ITEM(c,2,0), v4);
    element_set(ITEM(c,2,1), v5);
    element_clear(a0);
    element_clear(a1);
    element_clear(a2);
    element_clear(a3);
    element_clear(a4);
    element_clear(a5);
    element_clear(a6);
    element_clear(m0);
    element_clear(m1);
    element_clear(m2);
    element_clear(m3);
    element_clear(m4);
    element_clear(m5);
    element_clear(m6);
    element_clear(m7);
    element_clear(m8);
    element_clear(v0);
    element_clear(v1);
    element_clear(v2);
    element_clear(v3);
    element_clear(v4);
    element_clear(v5);
    element_clear(e1);
}

/* computing $c <- U^M, M=(3^{3m}-1)*(3^m+1)*(3^m+1-\mu*b*3^{(m+1)//2})$
 * This is the algorithm 8 in the paper above. */
static void algorithm8(element_t c, element_t u) {
    field_ptr f6 = FIELD(u), f = FIELD(ITEM(u,0,0));
    params *p = (params *) f->data;
    element_t v, w;
    element_init(v, f6);
    element_init(w, f6);
    algorithm6(v, u);
    algorithm7(v, v);
    element_set(w, v);
    int i;
    for (i = 0; i < (p->m + 1) / 2; i++)
        element_cubic(w, w);
    algorithm7(v, v);
    if (p->m % 12 == 1 || p->m % 12 == 11) { // w <= w^{-\mu*b}
        element_ptr e;
        e = ITEM(w,0,1);
        element_neg(e, e);
        e = ITEM(w,1,1);
        element_neg(e, e);
        e = ITEM(w,2,1);
        element_neg(e, e);
    }
    element_mul(c, v, w);
    element_clear(v);
    element_clear(w);
}

/* computing the Eta_T bilinear pairing $c <- Eta_T pairing(P,R)$ */
static void eta_T_pairing(element_ptr c, element_ptr P, element_ptr R, struct pairing_s *p) {
    UNUSED_VAR(p);
    if (DATA(P)->isinf || DATA(R)->isinf)
        element_set1(c);
    else {
        element_ptr x1 = DATA(P)->x, y1 = DATA(P)->y, x2 = DATA(R)->x, y2 =
                DATA(R)->y;
        if((PARAM(x1)->m - 1) / 2 % 2 == 0)
            algorithm5(c, x1, y1, x2, y2);
        else
            algorithm4(c, x1, y1, x2, y2);
        algorithm8(c, c);
    }
}

static void eta_T_3_clear(params *p) {
    mpz_clear(p->n);
    mpz_clear(p->n2);
    pbc_free(p);
}

static void GT_random(element_ptr e) {
    element_t a, b;
    element_init(a, e->field->pairing->G1);
    element_init(b, e->field->pairing->G1);
    element_random(a);
    element_random(b);
    element_pairing(e, a, b);
    element_clear(a);
    element_clear(b);
}

static void eta_T_3_pairing_clear(pairing_t pairing) {
    mpz_clear(pairing->r);
    field_clear(pairing->Zr);
    field_clear(pairing->GT);
    field_clear(pairing->G1);
    pbc_free(pairing->G1);
    pairing_data_ptr dp = pairing->data;
    field_clear(dp->gf3m);
    field_clear(dp->gf32m);
    field_clear(dp->gf36m);
    mpz_clear(dp->n2);
    pbc_free(dp);
}

static void eta_T_3_init_pairing(pairing_t pairing, params *p) {
    mpz_init(pairing->r);
    mpz_set(pairing->r, p->n);
    field_init_fp(pairing->Zr, pairing->r);

    pairing_data_ptr dp = pbc_malloc(sizeof(*dp));
    mpz_init(dp->n2);
    mpz_set(dp->n2, p->n2);
    field_init_gf3m(dp->gf3m, p->m, p->t);
    field_init_gf32m(dp->gf32m, dp->gf3m);
    field_init_gf33m(dp->gf36m, dp->gf32m);
    pairing_GT_init(pairing, dp->gf36m);
    pairing->GT->name = "eta_T_3 group of roots of 1";
    pairing->GT->random = GT_random;
    pairing->G2 = pairing->G1 = pbc_malloc(sizeof(field_t));
    field_init_eta_T_3(pairing->G1, dp->gf3m);
    pairing->G1->pairing = pairing;
    mpz_set(pairing->G1->order, p->n);
    mpz_set(pairing->GT->order, p->n);
    pairing->map = eta_T_pairing;
    pairing->data = dp;
    pairing->clear_func = eta_T_3_pairing_clear;
}

static void eta_T_3_out_str(FILE *stream, params *p) {
    param_out_type(stream, "i");
    param_out_int(stream, "m", p->m);
    param_out_int(stream, "t", p->t);
    param_out_mpz(stream, "n", p->n);
    param_out_mpz(stream, "n2", p->n2);
}

static void param_init(pbc_param_ptr p) {
    static pbc_param_interface_t interface = {{
      (void (*)(void *))eta_T_3_clear,
      (void (*)(pairing_t, void *))eta_T_3_init_pairing,
      (void (*)(FILE *, void *))eta_T_3_out_str,
    }};
    p->api = interface;
    params *param = p->data = pbc_malloc(sizeof(*param));
    mpz_init(param->n);
    mpz_init(param->n2);
}

int pbc_param_init_i(pbc_param_ptr p, struct symtab_s *tab) {
    param_init(p);
    params *param = p->data;
    int err = 0;
    err += lookup_int(&param->m, tab, "m");
    err += lookup_int(&param->t, tab, "t");
    err += lookup_mpz(param->n, tab, "n");
    err += lookup_mpz(param->n2, tab, "n2");
    return err;
}

void pbc_param_init_i_gen(pbc_param_ptr par, int group_size) {
    param_init(par);
    params *p = par->data;
    if (group_size <= 150) {
        p->m = 97;
        p->t = 12;
        mpz_set_str(p->n, "2726865189058261010774960798134976187171462721", 10);
        mpz_set_str(p->n2, "7", 10);
    } else if (group_size <= 206) {
        p->m = 199;
        p->t = 164;
        mpz_set_str(p->n, "167725321489096000055336949742738378351010268990525380470313869", 10);
        mpz_set_str(p->n2, "527874953560391326545598291952743", 10);
    } else if (group_size <= 259) {
        p->m = 235;
        p->t = 26;
        mpz_set_str(p->n, "1124316700897695330265827797088699345032488681307846555184025129863722718180241", 10);
        mpz_set_str(p->n2, "11819693021332914275777073321995059", 10);
    } else if (group_size <= 316) {
        p->m = 385;
        p->t = 22;
        mpz_set_str(p->n, "140884762419712839999909157778648717913595360839856026704744558309545986970238264714753014287541", 10);
        mpz_set_str(p->n2, "34899486997246711147841377458771182755186809219564106252058066150110543296498189654810187", 10);
    } else if (group_size <= 376) {
        p->m = 337;
        p->t = 30;
        mpz_set_str(p->n, "250796519030408069744426774377542635685621984993105288007781750196791322190409525696108840742205849171229571431053", 10);
        mpz_set_str(p->n2, "245777055088325363697128811262733732423405120899", 10);
    } else if (group_size <= 430) {
        p->m = 373;
        p->t = 198;
        mpz_set_str(p->n, "2840685307599487500956683789051368080919805957805957356540760731597378326586402072132959867084691357708217739285576524329854284197", 10);
        mpz_set_str(p->n2, "3256903458766749542151641063558247849550904613763", 10);
    } else if (group_size <= 484) {
        p->m = 395;
        p->t = 338;
        mpz_set_str(p->n, "80172097064154181257340545445945701478615643539554910656655431171167598268341527430200810544156625333601812351266052856520678455274751591367269291", 10);
        mpz_set_str(p->n2, "3621365590261279902324876775553649595261567", 10);
    } else if (group_size <= 552) {
        p->m = 433;
        p->t = 120;
        mpz_set_str(p->n, "15699907553631673835088720676147779193076555382157913339177784853763686462870506492752576492212322736133645158157557950634628006965882177348385366381692092784577773463", 10);
        mpz_set_str(p->n2, "24980791723059119877470531054938874784049", 10);
    } else if (group_size <= 644) {
        p->m = 467;
        p->t = 48;
        mpz_set_str(p->n, "108220469499363631995525712756135494735252733492048868417164002000654321383482753640072319529019505742300964525569770933946381504691909098938045089999753901375631613294579329433690943459352138231", 10);
        mpz_set_str(p->n2, "60438898450096967424971813347", 10);
    } else if (group_size <= 696) {
        p->m = 503;
        p->t = 104;
        mpz_set_str(p->n, "545523657676112447260904563578912738373307867219686215849632469801471112426878939776725222290437653718473962733760874627315930933126581248465899651120481066111839081575164964589811985885719017214938514563804313", 10);
        mpz_set_str(p->n2, "1799606423432800810122901025413", 10);
    } else if (group_size <= 803) {
        p->m = 509;
        p->t = 358;
        mpz_set_str(p->n, "102239946202586852409809887418093021457150612495255706614733003327526279081563687830782748305746187060264985869283524441819589592750998086186315250781067131293823177124077445718802216415539934838376431091001197641295264650596195201747790167311", 10);
        mpz_set_str(p->n2, "7", 10);
    } else if (group_size <= 892) {
        p->m = 617;
        p->t = 88;
        mpz_set_str(p->n, "57591959284219511220590893724691916802833742568034971006633345422620650391172287893878655658086794200963521584019889327992536532560877385225451713282279597074750857647455565899702728629166541223955196002755787520206774906606158388947359746178875040401304783332742806641", 10);
        mpz_set_str(p->n2, "42019638181715250622338241", 10);
    } else
        pbc_die("unsupported group size");
}

