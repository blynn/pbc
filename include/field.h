/*
 * Poorly named. Though this data structure can represent a field,
 * it is more general, as it is used to implement rings and groups as well.
 */
#ifndef FIELD_H
#define FIELD_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "random.h"

struct field_s;

struct element_s {
    struct field_s *field;
    void *data;
};
typedef struct element_s *element_ptr;
typedef struct element_s element_t[1];

struct field_s {
    void (*field_clear)(struct field_s *f);
    void (*init)(element_ptr);
    void (*clear)(element_ptr);
    void (*set_si)(element_ptr, signed long int);
    void (*set_mpz)(element_ptr, mpz_ptr);
    void (*set)(element_ptr, element_ptr);
    void (*set0)(element_ptr);
    void (*set1)(element_ptr);
    size_t (*out_str)(FILE *stream, int base, element_ptr);
    void (*add)(element_ptr, element_ptr, element_ptr);
    void (*sub)(element_ptr, element_ptr, element_ptr);
    void (*mul)(element_ptr, element_ptr, element_ptr);
    void (*mul_mpz)(element_ptr, element_ptr, mpz_ptr);
    void (*mul_si)(element_ptr, element_ptr, signed long int);
    void (*square)(element_ptr, element_ptr);
    void (*pow)(element_ptr, element_ptr, mpz_ptr);
    void (*invert)(element_ptr, element_ptr);
    void (*neg)(element_ptr, element_ptr);
    int (*cmp)(element_ptr, element_ptr);
    void (*random)(element_ptr);
    void (*from_hash)(element_ptr, int len, void *data);
    int (*is1)(element_ptr);
    int (*is0)(element_ptr);
    int (*sign)(element_ptr); //satisfies sign(x) = -sign(-x)
    int (*is_sqr)(element_ptr);
    void (*sqrt)(element_ptr, element_ptr);
    int (*to_bytes)(unsigned char *data, element_ptr);
    int (*from_bytes)(element_ptr, unsigned char *data);
    int (*length_in_bytes)(element_ptr);
    int fixed_length_in_bytes; //length of an element in bytes; -1 for variable
    mpz_t order; //-1 for infinite order
    element_ptr nqr; //nonquadratic residue
    void (*to_mpz)(mpz_ptr, element_ptr);
    void *data;
};
typedef struct field_s *field_ptr;
typedef struct field_s field_t[1];

typedef void (*fieldmap)(element_t dst, element_t src);

/*@manual internal
Initialize ''e'' to be an element of the group, ring or field ''f''
and set it to be the zero element.
*/
static inline void element_init(element_t e, field_ptr f)
{
    e->field = f;
    f->init(e);
}

/*@manual einit
Free the space occupied by ''e''. Call this when
the variable ''e'' is no longer needed.
*/
static inline void element_clear(element_t e)
{
    e->field->clear(e);
}

static inline size_t element_out_str(FILE *stream, int base, element_t e)
{
    return e->field->out_str(stream, base, e);
}

int element_fprintf(FILE *stream, const char *format, ...);
int element_printf(const char *format, ...);

/*@manual eassign
Set ''e'' to zero.
*/
static inline void element_set0(element_t e)
{
    e->field->set0(e);
}

/*@manual eassign
Set ''e'' to one.
*/
static inline void element_set1(element_t e)
{
    e->field->set1(e);
}

/*@manual eassign
Set ''e'' to ''i''.
*/
static inline void element_set_si(element_t e, signed long int i)
{
    e->field->set_si(e, i);
}

/*@manual eassign
Set ''e'' to ''z''.
*/
static inline void element_set_mpz(element_t e, mpz_t z)
{
    e->field->set_mpz(e, z);
}

/*@manual eassign
Set ''e'' to ''a''.
*/
static inline void element_set(element_t e, element_t a)
{
    e->field->set(e, a);
}

/*@manual earith
Set ''n'' to ''a'' + ''b''.
*/
static inline void element_add(element_t n, element_t a, element_t b)
{
    n->field->add(n, a, b);
}

/*@manual earith
Set ''n'' to ''a'' - ''b''.
*/
static inline void element_sub(element_t n, element_t a, element_t b)
{
    n->field->sub(n, a, b);
}

/*@manual earith
Set ''n'' to ''a'' times ''b''.
*/
static inline void element_mul(element_t n, element_t a, element_t b)
{
    n->field->mul(n, a, b);
}

/*@manual earith
*/
static inline void element_mul_mpz(element_t n, element_t a, mpz_t z)
{
    n->field->mul_mpz(n, a, z);
}

/*@manual earith
Set ''n'' to ''a'' times ''z'', that is ''a'' + ''a'' + ... + ''a''
where there are ''z'' ''a'''s.
*/
static inline void element_mul_si(element_t n, element_t a, signed long int z) 
{
    n->field->mul_si(n, a, z);
}

/*@manual earith
Set ''n'' to ''a'' times ''a''.
*/
static inline void element_square(element_t n, element_t a)
{
    n->field->square(n, a);
}

/*@manual earith
Set ''n'' to ''a'' raised to the power ''exp'', that is
''a'' times ''a'' times ... times ''a'' where there are ''exp'' ''a'''s
*/
static inline void element_pow(element_t x, element_t a, mpz_t exp)
{
    x->field->pow(x, a, exp);
}

static inline void element_pow_fp(element_t x, element_t a, element_t n)
    //n should be an element of F_p
    //(p should be the order of x)
{
    x->field->pow(x, a, n->data);
}

/*@manual earith
Set ''n'' to -''a''.
*/
static inline void element_neg(element_t n, element_t a)
{
    n->field->neg(n, a);
}

/*@manual earith
Set ''n'' to the inverse of ''a''.
*/
static inline void element_invert(element_t n, element_t a)
{
    n->field->invert(n, a);
}

static inline void element_random(element_t n)
{
    n->field->random(n);
}

/*@manual ecmp
Returns 0 if ''n'' is 1, nonzero otherwise.
*/
static inline int element_is1(element_t n)
{
    return n->field->is1(n);
}

/*@manual ecmp
Returns 0 if ''n'' is 0, nonzero otherwise.
*/
static inline int element_is0(element_t n)
{
    return n->field->is0(n);
}

/*@manual ecmp
Returns 0 if ''a'' and ''b are the same, nonzero otherwise.
*/
static inline int element_cmp(element_t a, element_t b)
{
    return a->field->cmp(a, b);
}

/*@manual ecmp
Returns 0 if ''a'' is a perfect square (quadratic residue), nonzero otherwise.
*/
static inline int element_is_sqr(element_t a)
{
    return a->field->is_sqr(a);
}

/*@manual ecmp
*/
static inline int element_sgn(element_t a)
{
    return a->field->sign(a);
}

/*@manual ecmp
If ''a'' is zero, returns 0. Otherwise:
For ''a'' in a field GF(p), returns -1 if ''a'' &lt; p, 1 otherwise.
For ''a'' in a polynomial ring, returns <function>element_sgn</function>
called on the coefficient of the lowest degree term.
*/
static inline int element_sign(element_t a)
{
    return a->field->sign(a);
}

static inline void element_sqrt(element_t a, element_t b)
{
    a->field->sqrt(a, b);
}

static inline void element_from_hash(element_t a, int len, void *data)
{
    a->field->from_hash(a, len, data);
}

static inline int element_to_bytes(unsigned char *data, element_t a)
{
    return a->field->to_bytes(data, a);
}

static inline int element_from_bytes(element_t a, unsigned char *data)
{
    return a->field->from_bytes(a, data);
}

/*@manual econvert
Converts ''e'' to a GMP integer ''z''
if such an operation makes sense
*/
static inline void element_to_mpz(mpz_t z, element_t e)
{
    e->field->to_mpz(z, e);
}

static inline int element_length_in_bytes(element_t a)
{
    if (a->field->fixed_length_in_bytes < 0) {
	return a->field->length_in_bytes(a);
    } else {
	return a->field->fixed_length_in_bytes;
    }
}

void element_pow2(element_t x, element_t a1, mpz_t n1,
                                 element_t a2, mpz_t n2);
void element_pow3(element_t x, element_t a1, mpz_t n1,
                                 element_t a2, mpz_t n2,
                                 element_t a3, mpz_t n3);

static inline void field_clear(field_ptr f)
{
    if (f->nqr) {
	element_clear(f->nqr);
	free(f->nqr);
    }
    mpz_init(f->order);
    f->field_clear(f);
}

element_ptr field_get_nqr(field_ptr f);

void field_init(field_ptr f);

static inline int mpz_is0(mpz_t z)
{
    return !mpz_sgn(z);
    //return !mpz_cmp_ui(z, 0);
}

void field_init_naive_fp(field_ptr f, mpz_t prime);
void field_init_slow_fp(field_ptr f, mpz_t prime);
void field_init_fast_fp(field_ptr f, mpz_t prime);

static inline void field_init_fp(field_ptr f, mpz_t prime) {
    field_init_naive_fp(f, prime);
}

//The following only work for groups of points on elliptic curves:
int element_to_bytes_x_only(unsigned char *data, element_t e);
int element_from_bytes_x_only(element_t e, unsigned char *data);
int element_length_in_bytes_x_only(element_t e);
int element_to_bytes_compressed(unsigned char *data, element_t e);
int element_from_bytes_compressed(element_t e, unsigned char *data);
int element_length_in_bytes_compressed(element_t e);

#endif //FIELD_H
