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
    void (*to_mpz)(mpz_ptr, element_ptr); //does this make sense for all fields?
    void *data;
};
typedef struct field_s *field_ptr;
typedef struct field_s field_t[1];

typedef void (*fieldmap)(element_t dst, element_t src);

static inline void element_init(element_ptr e, field_ptr f)
{
    e->field = f;
    f->init(e);
}

static inline void element_clear(element_ptr e)
{
    e->field->clear(e);
    //e->field = NULL;
}

static inline size_t element_out_str(FILE *stream, int base, element_ptr e)
{
    return e->field->out_str(stream, base, e);
}

int element_fprintf(FILE *stream, const char *format, ...);
int element_printf(const char *format, ...);

static inline void element_set_si(element_ptr e, signed long int i)
{
    e->field->set_si(e, i);
}

static inline void element_set_mpz(element_ptr e, mpz_ptr z)
{
    e->field->set_mpz(e, z);
}

static inline void element_set0(element_ptr e)
{
    e->field->set0(e);
}

static inline void element_set1(element_ptr e)
{
    e->field->set1(e);
}

static inline void element_set(element_ptr x, element_ptr a)
{
    x->field->set(x, a);
}

static inline void element_add(element_ptr n, element_ptr a, element_ptr b)
{
    n->field->add(n, a, b);
}

static inline void element_sub(element_ptr n, element_ptr a, element_ptr b)
{
    n->field->sub(n, a, b);
}

static inline void element_mul(element_ptr n, element_ptr a, element_ptr b)
{
    n->field->mul(n, a, b);
}

static inline void element_mul_mpz(element_ptr n, element_ptr a, mpz_ptr z)
{
    n->field->mul_mpz(n, a, z);
}

static inline void element_mul_si(element_ptr n, element_ptr a, signed long int z) 
{
    n->field->mul_si(n, a, z);
}

static inline void element_square(element_ptr n, element_ptr a)
{
    n->field->square(n, a);
}

static inline void element_pow(element_ptr x, element_ptr a, mpz_ptr n)
{
    x->field->pow(x, a, n);
}

static inline void element_pow_fp(element_ptr x, element_ptr a, element_ptr n)
    //n should be an element of F_p
    //(p should be the order of x)
{
    x->field->pow(x, a, n->data);
}

static inline void element_neg(element_ptr n, element_ptr a)
{
    n->field->neg(n, a);
}

static inline void element_invert(element_ptr n, element_ptr a)
{
    n->field->invert(n, a);
}

static inline void element_random(element_ptr n)
{
    n->field->random(n);
}

static inline int element_is1(element_ptr n)
{
    return n->field->is1(n);
}

static inline int element_is0(element_ptr n)
{
    return n->field->is0(n);
}

static inline int element_cmp(element_ptr a, element_ptr b)
{
    return a->field->cmp(a, b);
}

static inline int element_is_sqr(element_ptr a)
{
    return a->field->is_sqr(a);
}

static inline int element_sign(element_ptr a)
{
    return a->field->sign(a);
}

static inline void element_sqrt(element_ptr a, element_ptr b)
{
    a->field->sqrt(a, b);
}

static inline void element_from_hash(element_ptr a, int len, void *data)
{
    a->field->from_hash(a, len, data);
}

static inline int element_to_bytes(unsigned char *data, element_ptr a)
{
    return a->field->to_bytes(data, a);
}

static inline int element_from_bytes(element_ptr a, unsigned char *data)
{
    return a->field->from_bytes(a, data);
}

static inline void element_to_mpz(mpz_ptr z, element_ptr a)
{
    a->field->to_mpz(z, a);
}

static inline int element_length_in_bytes(element_ptr a)
{
    if (a->field->fixed_length_in_bytes < 0) {
	return a->field->length_in_bytes(a);
    } else {
	return a->field->fixed_length_in_bytes;
    }
}

void element_pow2(element_ptr x, element_ptr a1, mpz_ptr n1,
                                 element_ptr a2, mpz_ptr n2);
void element_pow3(element_ptr x, element_ptr a1, mpz_ptr n1,
                                 element_ptr a2, mpz_ptr n2,
                                 element_ptr a3, mpz_ptr n3);

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

static inline int mpz_is0(mpz_ptr z)
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
int element_to_bytes_x_only(unsigned char *data, element_ptr e);
int element_from_bytes_x_only(element_ptr e, unsigned char *data);
int element_length_in_bytes_x_only(element_ptr e);
int element_to_bytes_compressed(unsigned char *data, element_ptr e);
int element_from_bytes_compressed(element_ptr e, unsigned char *data);
int element_length_in_bytes_compressed(element_ptr e);

#endif //FIELD_H
