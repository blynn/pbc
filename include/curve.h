#ifndef CURVE_H
#define CURVE_H

#include "poly.h"

struct curve_s;

struct point_s {
    int inf_flag;
    element_t x;
    element_t y;
    struct curve_s *curve;
};
typedef struct point_s *point_ptr;
typedef struct point_s point_t[1];

struct curve_s {
    field_ptr field;
    void (*random)(point_ptr);
    void (*from_x)(point_ptr, element_ptr);
    void (*from_hash)(point_ptr, int, void *);
    void (*neg)(point_ptr, point_ptr);
    void (*add)(point_ptr, point_ptr, point_ptr);
    void (*double_nocheck)(point_ptr, point_ptr);
    void (*doublefn)(point_ptr, point_ptr);
    void (*mul)(point_ptr, mpz_ptr, point_ptr);
    void (*curve_clear)(struct curve_s *curve);
    void *data;
};
typedef struct curve_s *curve_ptr;
typedef struct curve_s curve_t[1];

static inline void curve_clear(curve_ptr c)
{
    c->curve_clear(c);
}

static inline void point_set_inf(point_ptr p)
{
    p->inf_flag = 1;
}

static inline void point_init(point_ptr p, curve_ptr c)
{
    field_ptr f = c->field;
    element_init(p->x, f);
    element_init(p->y, f);
    p->curve = c;
    point_set_inf(p);
}

static inline void point_clear(point_ptr p)
{
    //p->curve = NULL;
    element_clear(p->x);
    element_clear(p->y);
}

static inline void point_random(point_ptr p)
{
    p->curve->random(p);
}

static inline void point_from_hash(point_ptr p, int len, void *data)
{
    p->curve->from_hash(p, len, data);
}

static inline void point_neg(point_ptr r, point_ptr p)
{
    p->curve->neg(r, p);
}

static inline void point_add(point_ptr r, point_ptr p, point_ptr q)
{
    r->curve->add(r, p, q);
}

static inline void point_double(point_ptr r, point_ptr p)
{
    r->curve->doublefn(r, p);
}

static inline void point_set(point_ptr p, point_ptr q)
{
    if (q->inf_flag) {
	p->inf_flag = 1;
	return;
    }
    p->inf_flag = 0;
    element_set(p->x, q->x);
    element_set(p->y, q->y);
}

size_t point_out_str(FILE *stream, int base, point_ptr p);

static inline void point_mul(point_ptr r, mpz_ptr n, point_ptr p)
{
    r->curve->mul(r, n, p);
}

static inline int point_is_inf(point_ptr p)
{
    return p->inf_flag;
}

static inline void point_from_x(point_ptr p, element_ptr x)
{
    p->curve->from_x(p, x);
}

void point_map(point_t R, fieldmap map, point_t P);

//I define "common" curves (abbr. "cc") to be those
//of the form y^2 = x^3 + ax + b defined over a field with
//characteristic > 3

struct common_curve_s {
    element_t a, b;
};
typedef struct common_curve_s common_curve_t[1];
typedef struct common_curve_s *common_curve_ptr;

//for shoehorning the group of points on a curve into the `field' data type
struct curve_group_data_s {
    curve_ptr curve;
    mpz_t cofac;
};
typedef struct curve_group_data_s curve_group_data_t[1];
typedef struct curve_group_data_s *curve_group_data_ptr;

/* called in {e,f}_param.c */
void curve_init_b(curve_ptr c, element_ptr b);

/* these are called in mnt.c */
void curve_init_cc_j(curve_ptr c, element_ptr j);
    //assumes j != 0, 1728
void twist_curve(curve_ptr c);
    //WARNING: existing points will no longer lie on c
    //as this modifies c in place
void compute_trace_n(mpz_t res, mpz_t q, mpz_t trace, int n);
    //compute trace of Frobenius at q^n given trace at q
    //see p.105 of Blake, Seroussi and Smart

/* these are called in pairing.c */
void cc_frobenius(point_ptr r, point_ptr p, mpz_ptr q);
void curve_init_cc_ab(curve_ptr c, element_ptr a, element_ptr b);
void cc_init_map_curve(curve_ptr cnew, curve_ptr c,
                       field_ptr dstfield, fieldmap map);
void field_init_curve_group(field_t f, curve_t c, mpz_t cofac);
void curve_init_singular_with_node(curve_ptr c, field_t field);

#endif //CURVE_H
