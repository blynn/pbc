/*
 * Though this data structure does represent fields,
 * we also use it for rings and groups as well.
 */
//requires
// * stdio.h
// * gmp.h
#ifndef __PBC_FIELD_H__
#define __PBC_FIELD_H__

struct field_s;

struct element_s {
    struct field_s *field;
    void *data;
};
typedef struct element_s *element_ptr;
typedef struct element_s element_t[1];

struct element_pp_s {
    struct field_s *field;
    void *data;
};
typedef struct element_pp_s element_pp_t[1];
typedef struct element_pp_s *element_pp_ptr;

#include "pbc_assert.h"
void pbc_assert(int expr, char *msg, const char *func);
void pbc_assert_match2(element_ptr a, element_ptr b, const char *func);
void pbc_assert_match3(element_ptr a, element_ptr b, element_ptr c, const char *func);

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
    void (*add_ui)(element_ptr, element_ptr, unsigned long int);
    void (*sub)(element_ptr, element_ptr, element_ptr);
    void (*mul)(element_ptr, element_ptr, element_ptr);
    void (*mul_mpz)(element_ptr, element_ptr, mpz_ptr);
    void (*mul_si)(element_ptr, element_ptr, signed long int);
    void (*div)(element_ptr, element_ptr, element_ptr);
    void (*square)(element_ptr, element_ptr);
    int (*is_sqr)(element_ptr);
    void (*sqrt)(element_ptr, element_ptr);
    void (*doub)(element_ptr, element_ptr); //can't call it "double"!
    void (*halve)(element_ptr, element_ptr);
    void (*pow_mpz)(element_ptr, element_ptr, mpz_ptr);
    void (*invert)(element_ptr, element_ptr);
    void (*neg)(element_ptr, element_ptr);
    void (*random)(element_ptr);
    void (*from_hash)(element_ptr, void *data, int len);
    int (*is1)(element_ptr);
    int (*is0)(element_ptr);
    int (*sign)(element_ptr); //satisfies sign(x) = -sign(-x)
    int (*cmp)(element_ptr, element_ptr);
    int (*to_bytes)(unsigned char *data, element_ptr);
    int (*from_bytes)(element_ptr, unsigned char *data);
    int (*length_in_bytes)(element_ptr);
    int fixed_length_in_bytes; //length of an element in bytes; -1 for variable
    mpz_t order; //0 for infinite order
    element_ptr nqr; //nonquadratic residue
    void (*to_mpz)(mpz_ptr, element_ptr);
    void (*out_info)(FILE *, struct field_s *);
    void (*pp_init)(element_pp_t p, element_t in);
    void (*pp_clear)(element_pp_t p);
    void (*pp_pow)(element_t out, mpz_ptr power, element_pp_t p);
    int (*snprint)(char *s, size_t n, element_ptr e);
    void *data;
};
typedef struct field_s *field_ptr;
typedef struct field_s field_t[1];

typedef void (*fieldmap)(element_t dst, element_t src);

/*@manual internal
Initialize ''e'' to be an element of the algebraic structure ''f''
and set it to be the zero element.
*/
static inline void element_init(element_t e, field_ptr f)
{
    e->field = f;
    f->init(e);
}

/*@manual einit
Initialize ''e'' to be an element of the algebraic structure that ''e2''
lies in.
*/
static inline void element_init_same_as(element_t e, element_t e2)
{
    element_init(e, e2->field);
}

/*@manual einit
Free the space occupied by ''e''. Call this when
the variable ''e'' is no longer needed.
*/
static inline void element_clear(element_t e)
{
    e->field->clear(e);
}

/*@manual eio
Output ''e'' on ''stream'' in base ''base''. The base must be between
2 and 36.
*/
static inline size_t element_out_str(FILE *stream, int base, element_t e)
{
    return e->field->out_str(stream, base, e);
}

/*@manual eio
*/
int element_printf(const char *format, ...);

/*@manual eio
*/
int element_fprintf(FILE *stream, const char *format, ...);

/*@manual eio
Same as printf family
except also has the 'B' conversion specifier for types
of <function>element_t</function>, and 'Y', 'Z' conversion specifiers for
<type>mpz_t</type>. For example if ''e'' is of type
<type>element_t</type> then
<screen>
element_printf("%B\n", e);
</screen>
will print the value of ''e'' in a human-readable form on standard output.
*/
int element_snprintf (char *buf, size_t size, const char *fmt, ...);

/*@manual eio
Convert an element to a human-friendly string.
Behaves as <function>snprintf</function> but only on one element at a time.
*/
static inline int element_snprint(char *s, size_t n, element_t e)
{
    return e->field->snprint(s, n, e);
}

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
    PBC_ASSERT_MATCH2(e, a);
    e->field->set(e, a);
}

static inline void element_add_ui(element_t n, element_t a, unsigned long int b)
{
    n->field->add_ui(n, a, b);
}

/*@manual econvert
Converts ''e'' to a GMP integer ''z''
if such an operation makes sense
*/
static inline void element_to_mpz(mpz_t z, element_t e)
{
    e->field->to_mpz(z, e);
}

/*@manual econvert
Generate an element ''e'' deterministically from
the ''len'' bytes stored in the buffer ''data''.
*/
static inline void element_from_hash(element_t e, void *data, int len)
{
    e->field->from_hash(e, data, len);
}
/*@manual earith
Set ''n'' to ''a'' + ''b''.
*/
static inline void element_add(element_t n, element_t a, element_t b)
{
    PBC_ASSERT_MATCH3(n, a, b);
    n->field->add(n, a, b);
}

/*@manual earith
Set ''n'' to ''a'' - ''b''.
*/
static inline void element_sub(element_t n, element_t a, element_t b)
{
    PBC_ASSERT_MATCH3(n, a, b);
    n->field->sub(n, a, b);
}

/*@manual earith
Set ''n'' to ''a'' times ''b''.
*/
static inline void element_mul(element_t n, element_t a, element_t b)
{
    PBC_ASSERT_MATCH3(n, a, b);
    n->field->mul(n, a, b);
}

/*@manual earith
*/
static inline void element_mul_mpz(element_t n, element_t a, mpz_t z)
{
    PBC_ASSERT_MATCH2(n, a);
    n->field->mul_mpz(n, a, z);
}

/*@manual earith
''z'' must be a mpz_t/signed long int.
Set ''n'' to ''a'' times ''z'', that is ''a'' + ''a'' + ... + ''a''
where there are ''z'' ''a'''s.
*/
static inline void element_mul_si(element_t n, element_t a, signed long int z) 
{
    PBC_ASSERT_MATCH2(n, a);
    n->field->mul_si(n, a, z);
}

/*@manual earith
''z'' must be an element of a integer mod ring (i.e. Z_n for some n).
Set ''c'' to ''a'' times ''z'', that is ''a'' + ''a'' + ... + ''a''
where there are ''z'' ''a'''s.
*/
static inline void element_mul_zn(element_t c, element_t a, element_t z)
{
    mpz_t z0;
    PBC_ASSERT_MATCH2(c, a);
    //TODO: check z->field is Zn
    mpz_init(z0);
    element_to_mpz(z0, z);
    element_mul_mpz(c, a, z0);
    mpz_clear(z0);
}

/*@manual earith
Set ''n'' to ''a'' / ''b''.
*/
static inline void element_div(element_t n, element_t a, element_t b)
{
    PBC_ASSERT_MATCH3(n, a, b);
    n->field->div(n, a, b);
}

/*@manual earith
Set ''n'' to ''a'' + ''a''.
*/
static inline void element_double(element_t n, element_t a)
{
    PBC_ASSERT_MATCH2(n, a);
    n->field->doub(n, a);
}

/*@manual earith
Set ''n'' to ''a/2''
*/
static inline void element_halve(element_t n, element_t a)
{
    PBC_ASSERT_MATCH2(n, a);
    n->field->halve(n, a);
}

/*@manual earith
Set ''n'' to ''a'' times ''a''.
*/
static inline void element_square(element_t n, element_t a)
{
    PBC_ASSERT_MATCH2(n, a);
    n->field->square(n, a);
}

/*@manual epow
Set ''x'' to ''a'' raised to the power ''n'', that is
''a'' times ''a'' times ... times ''a'' where there are ''n'' ''a'''s
*/
static inline void element_pow_mpz(element_t x, element_t a, mpz_t n)
{
    PBC_ASSERT_MATCH2(x, a);
    x->field->pow_mpz(x, a, n);
}

/*@manual epow
Set ''x'' to ''a'' raised to the power ''n'', where ''n'' is
an element of a ring Z_n for some n (typically the order
of the algebraic structure ''x'' lies in).
*/
static inline void element_pow_zn(element_t x, element_t a, element_t n)
{
    mpz_t z;
    PBC_ASSERT_MATCH2(x, a);
    mpz_init(z);
    element_to_mpz(z, n);
    x->field->pow_mpz(x, a, z);
    mpz_clear(z);
}

/*@manual earith
Set ''n'' to -''a''.
*/
static inline void element_neg(element_t n, element_t a)
{
    PBC_ASSERT_MATCH2(n, a);
    n->field->neg(n, a);
}

/*@manual earith
Set ''n'' to the inverse of ''a''.
*/
static inline void element_invert(element_t n, element_t a)
{
    PBC_ASSERT_MATCH2(n, a);
    n->field->invert(n, a);
}

/*@manual erandom
If the ''e'' lies in a finite algebraic structure,
this function assigns a uniformly random element to ''e''.
*/
static inline void element_random(element_t e)
{
    e->field->random(e);
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
Returns 0 if ''a'' and ''b'' are the same, nonzero otherwise.
*/
static inline int element_cmp(element_t a, element_t b)
{
    PBC_ASSERT_MATCH2(a, b);
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
If ''a'' is zero, returns 0. For nozero ''a'' the behaviour depends on
the algebraic structure, but has the property that
element_sgn(''a'') = -element_sgn(-''a'')
and
element_sgn(''a'') = 0 implies ''a'' = 0 with overwhelming probability.
Not implemented on elliptic curve groups yet.
*/
static inline int element_sign(element_t a)
{
    return a->field->sign(a);
}

static inline void element_sqrt(element_t a, element_t b)
{
    PBC_ASSERT_MATCH2(a, b);
    a->field->sqrt(a, b);
}

/*@manual etrade
Returns the length in bytes the element ''e'' will take to represent
*/
static inline int element_length_in_bytes(element_t e)
{
    if (e->field->fixed_length_in_bytes < 0) {
	return e->field->length_in_bytes(e);
    } else {
	return e->field->fixed_length_in_bytes;
    }
}

/*@manual etrade
Converts ''e'' to byte, writing the result in the buffer ''data''.
The number of bytes it will write can be determined from calling
<function>element_length_in_bytes()</function>.
Returns number of bytes written.
*/
static inline int element_to_bytes(unsigned char *data, element_t e)
{
    return e->field->to_bytes(data, e);
}

/*@manual etrade
Reads ''e'' from the buffer ''data'', and returns
the number of bytes read.
*/
static inline int element_from_bytes(element_t e, unsigned char *data)
{
    return e->field->from_bytes(e, data);
}

/*@manual epow
Sets ''x'' = ''a1''^''n1'' times ''a2''^''n2'', and is generally faster than
performing two separate exponentiations.
*/
void element_pow2_mpz(element_t x, element_t a1, mpz_t n1,
                                 element_t a2, mpz_t n2);
/*@manual epow
Also sets ''x'' = ''a1''^''n1'' times ''a2''^''n2'',
but ''n1'', ''n2'' must be elements of a ring Z_n for some integer n.
*/
static inline void element_pow2_zn(element_t x, element_t a1, element_t n1,
                                 element_t a2, element_t n2)
{
    mpz_t z1, z2;
    mpz_init(z1);
    mpz_init(z2);
    element_to_mpz(z1, n1);
    element_to_mpz(z2, n2);
    element_pow2_mpz(x, a1, z1, a2, z2);
    mpz_clear(z1);
    mpz_clear(z2);
}

/*@manual epow
Sets ''x'' = ''a1''^''n1'' times ''a2^n2'' times ''a3''^''n3'',
and is generally faster than
performing three separate exponentiations.
*/
void element_pow3_mpz(element_t x, element_t a1, mpz_t n1,
                                 element_t a2, mpz_t n2,
                                 element_t a3, mpz_t n3);

/*@manual epow
Also sets ''x'' = ''a1''^''n1'' times ''a2^n2'' times ''a3''^''n3'',
but ''n1'', ''n2'', ''n3'' must be elements of a ring Z_n for some integer n.
*/
static inline void element_pow3_zn(element_t x, element_t a1, element_t n1,
                                 element_t a2, element_t n2,
                                 element_t a3, element_t n3)
{
    mpz_t z1, z2, z3;
    mpz_init(z1);
    mpz_init(z2);
    mpz_init(z3);
    element_to_mpz(z1, n1);
    element_to_mpz(z2, n2);
    element_to_mpz(z3, n3);
    element_pow3_mpz(x, a1, z1, a2, z2, a3, z3);
    mpz_clear(z1);
    mpz_clear(z2);
    mpz_clear(z3);
}

void field_clear(field_ptr f);

element_ptr field_get_nqr(field_ptr f);
void field_set_nqr(field_ptr f, element_t nqr);
void field_gen_nqr(field_ptr f);

void field_init(field_ptr f);

static inline int mpz_is0(mpz_t z)
{
    return !mpz_sgn(z);
    //return !mpz_cmp_ui(z, 0);
}

/*@manual etrade
Assumes ''e'' is a point on an elliptic curve.
Writes the x-coordinate of ''e'' to the buffer ''data''
*/
int element_to_bytes_x_only(unsigned char *data, element_t e);
/*@manual etrade
Assumes ''e'' is a point on an elliptic curve.
Sets ''e'' to a point with
x-coordinate represented by the buffer ''data''. This is not unique.
For each ''x''-coordinate, there exist two different points, at least
for the elliptic curves in PBC. (They are inverses of each other.)
*/
int element_from_bytes_x_only(element_t e, unsigned char *data);
/*@manual etrade
Assumes ''e'' is a point on an elliptic curve.
Returns the length in bytes needed to hold the x-coordinate of ''e''.
*/
int element_length_in_bytes_x_only(element_t e);

/*@manual etrade
If possible, outputs a compressed form of the element ''e'' to
the buffer of bytes ''data''.
Currently only implemented for points on an elliptic curve.
*/
int element_to_bytes_compressed(unsigned char *data, element_t e);

/*@manual etrade
Sets element ''e'' to the element in compressed form in the buffer of bytes
''data''.
Currently only implemented for points on an elliptic curve.
*/
int element_from_bytes_compressed(element_t e, unsigned char *data);

/*@manual etrade
Returns the number of bytes needed to hold ''e'' in compressed form.
Currently only implemented for points on an elliptic curve.
*/
int element_length_in_bytes_compressed(element_t e);

void field_out_info(FILE *out, field_ptr f);

/*@manual epow
Prepare to exponentiate an element ''in'', and store preprocessing information
in ''p''.
*/
static inline void element_pp_init(element_pp_t p, element_t in) {
    p->field = in->field;
    in->field->pp_init(p, in);
}

/*@manual epow
Clear ''p''. Should be called after ''p'' is no longer needed.
*/
static inline void element_pp_clear(element_pp_t p)
{
    p->field->pp_clear(p);
}

/*@manual epow
Raise ''in'' to ''power'' and store the result in ''out'', where ''in''
is a previously preprocessed element, that is, the second argument
passed to a previous <function>element_pp_init</function> call.
*/
static inline void element_pp_pow(element_t out, mpz_ptr power, element_pp_t p)
{
    p->field->pp_pow(out, power, p);
}

void pbc_mpz_out_raw_n(unsigned char *data, int n, mpz_t z);
void pbc_mpz_from_hash(mpz_t z, mpz_t limit,
	unsigned char *data, unsigned int len);

void brute_force_dlog(element_t x, element_t g, element_t h);
void pollard_rho(element_t x, element_t g, element_t h);
void index_calculus_dlog(mpz_t x, mpz_t g, mpz_t h, mpz_t q);

#endif //__PBC_FIELD_H__
