/*
 * field_t: represents fields, rings and groups.
 * element_t: represents an element of a field_t.
 */

// Requires:
// * stdarg.h
// * stdio.h
// * gmp.h
// * utils.h
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

void pbc_assert(int expr, char *msg, const char *func);
void pbc_assert_match2(element_ptr a, element_ptr b, const char *func);
void pbc_assert_match3(element_ptr a, element_ptr b, element_ptr c,
                       const char *func);

struct multiz_s;
typedef struct multiz_s *multiz;

struct pairing_s;
struct field_s {
  void (*field_clear)(struct field_s *f);
  void (*init)(element_ptr);
  void (*clear)(element_ptr);

  void (*set_mpz)(element_ptr, mpz_ptr);
  void (*set_multiz)(element_ptr, multiz);
  void (*set)(element_ptr, element_ptr);
  void (*set0)(element_ptr);
  void (*set1)(element_ptr);
  int (*set_str)(element_ptr e, const char *s, int base);
  size_t(*out_str)(FILE *stream, int base, element_ptr);
  void (*add)(element_ptr, element_ptr, element_ptr);
  void (*sub)(element_ptr, element_ptr, element_ptr);
  void (*mul)(element_ptr, element_ptr, element_ptr);

  int (*is_sqr)(element_ptr);
  void (*sqrt)(element_ptr, element_ptr);

  // Defaults exist for these functions.
  int (*item_count)(element_ptr);
  element_ptr (*item)(element_ptr, int);
  element_ptr (*get_x)(element_ptr);
  element_ptr (*get_y)(element_ptr);
  void (*set_si)(element_ptr, signed long int);
  void (*add_ui)(element_ptr, element_ptr, unsigned long int);
  void (*mul_mpz)(element_ptr, element_ptr, mpz_ptr);
  void (*mul_si)(element_ptr, element_ptr, signed long int);
  void (*div)(element_ptr, element_ptr, element_ptr);
  void (*doub)(element_ptr, element_ptr);  // Can't call it "double"!
  void (*multi_doub)(element_ptr*, element_ptr*, int n);
  void (*multi_add)(element_ptr*, element_ptr*, element_ptr*, int n);
  void (*halve)(element_ptr, element_ptr);
  void (*square)(element_ptr, element_ptr);

  void (*cubic) (element_ptr, element_ptr);
  void (*pow_mpz)(element_ptr, element_ptr, mpz_ptr);
  void (*invert)(element_ptr, element_ptr);
  void (*neg)(element_ptr, element_ptr);
  void (*random)(element_ptr);
  void (*from_hash)(element_ptr, void *data, int len);
  int (*is1)(element_ptr);
  int (*is0)(element_ptr);
  int (*sign)(element_ptr);  // satisfies sign(x) = -sign(-x)
  int (*cmp)(element_ptr, element_ptr);
  int (*to_bytes)(unsigned char *data, element_ptr);
  int (*from_bytes)(element_ptr, unsigned char *data);
  int (*length_in_bytes)(element_ptr);
  int fixed_length_in_bytes;  // length of an element in bytes; -1 for variable
  int (*snprint)(char *s, size_t n, element_ptr e);
  void (*to_mpz)(mpz_ptr, element_ptr);
  void (*out_info)(FILE *, struct field_s *);
  void (*pp_init)(element_pp_t p, element_t in);
  void (*pp_clear)(element_pp_t p);
  void (*pp_pow)(element_t out, mpz_ptr power, element_pp_t p);

  struct pairing_s *pairing;

  mpz_t order;                // 0 for infinite order
  element_ptr nqr;            // nonquadratic residue

  char *name;
  void *data;
};
typedef struct field_s *field_ptr;
typedef struct field_s field_t[1];

typedef void (*fieldmap) (element_t dst, element_t src);

void field_out_info(FILE* out, field_ptr f);

/*@manual internal
Initialize 'e' to be an element of the algebraic structure 'f'
and set it to be the zero element.
*/
static inline void element_init(element_t e, field_ptr f) {
  e->field = f;
  f->init(e);
}

element_ptr element_new(field_ptr f);
void element_free(element_ptr e);

/*@manual einit
Initialize 'e' to be an element of the algebraic structure that 'e2'
lies in.
*/
static inline void element_init_same_as(element_t e, element_t e2) {
  element_init(e, e2->field);
}

/*@manual einit
Free the space occupied by 'e'. Call this when
the variable 'e' is no longer needed.
*/
static inline void element_clear(element_t e) {
  e->field->clear(e);
}

/*@manual eio
Output 'e' on 'stream' in base 'base'. The base must be between
2 and 36.
*/
static inline size_t element_out_str(FILE * stream, int base, element_t e) {
  return e->field->out_str(stream, base, e);
}

/*@manual eio
*/
int element_printf(const char *format, ...);

/*@manual eio
*/
int element_fprintf(FILE * stream, const char *format, ...);

/*@manual eio
*/
int element_snprintf(char *buf, size_t size, const char *fmt, ...);

/*@manual eio
Same as printf family
except also has the 'B' conversion specifier for types
of *element_t*, and 'Y', 'Z' conversion specifiers for
+mpz_t+. For example if 'e' is of type
+element_t+ then

  element_printf("%B\n", e);

will print the value of 'e' in a human-readable form on standard output.
*/
int element_vsnprintf(char *buf, size_t size, const char *fmt, va_list ap);

/*@manual eio
Convert an element to a human-friendly string.
Behaves as *snprintf* but only on one element at a time.
*/
static inline int element_snprint(char *s, size_t n, element_t e) {
  return e->field->snprint(s, n, e);
}

static inline void element_set_multiz(element_t e, multiz m) {
  e->field->set_multiz(e, m);
}

/*@manual eio
Set the element 'e' from 's', a null-terminated C string in base 'base'.
Whitespace is ignored. Points have the form "['x,y']" or "'O'",
while polynomials have the form "['a0,...,an']".
Returns number of characters read (unlike GMP's mpz_set_str).
A return code of zero means PBC could not find a well-formed string
describing an element.
*/
static inline int element_set_str(element_t e, const char *s, int base) {
  return e->field->set_str(e, s, base);
}

/*@manual eassign
Set 'e' to zero.
*/
static inline void element_set0(element_t e) {
  e->field->set0(e);
}

/*@manual eassign
Set 'e' to one.
*/
static inline void element_set1(element_t e) {
  e->field->set1(e);
}

/*@manual eassign
Set 'e' to 'i'.
*/
static inline void element_set_si(element_t e, signed long int i) {
  e->field->set_si(e, i);
}

/*@manual eassign
Set 'e' to 'z'.
*/
static inline void element_set_mpz(element_t e, mpz_t z) {
  e->field->set_mpz(e, z);
}

/*@manual eassign
Set 'e' to 'a'.
*/
static inline void element_set(element_t e, element_t a) {
  PBC_ASSERT_MATCH2(e, a);
  e->field->set(e, a);
}

static inline void element_add_ui(element_t n, element_t a,
                                  unsigned long int b) {
  n->field->add_ui(n, a, b);
}

/*@manual econvert
Converts 'e' to a GMP integer 'z'
if such an operation makes sense
*/
static inline void element_to_mpz(mpz_t z, element_t e) {
  e->field->to_mpz(z, e);
}

static inline long element_to_si(element_t e) {
  mpz_t z;
  mpz_init(z);
  e->field->to_mpz(z, e);
  long res = mpz_get_si(z);
  mpz_clear(z);
  return res;
}

/*@manual econvert
Generate an element 'e' deterministically from
the 'len' bytes stored in the buffer 'data'.
*/
static inline void element_from_hash(element_t e, void *data, int len) {
  e->field->from_hash(e, data, len);
}

/*@manual earith
Set 'n' to 'a' + 'b'.
*/
static inline void element_add(element_t n, element_t a, element_t b) {
  PBC_ASSERT_MATCH3(n, a, b);
  n->field->add(n, a, b);
}

/*@manual earith
Set 'n' to 'a' - 'b'.
*/
static inline void element_sub(element_t n, element_t a, element_t b) {
  PBC_ASSERT_MATCH3(n, a, b);
  n->field->sub(n, a, b);
}

/*@manual earith
Set 'n' = 'a' 'b'.
*/
static inline void element_mul(element_t n, element_t a, element_t b) {
  PBC_ASSERT_MATCH3(n, a, b);
  n->field->mul(n, a, b);
}

static inline void element_cubic(element_t n, element_t a) {
  PBC_ASSERT_MATCH2(n, a);
  n->field->cubic(n, a);
}

/*@manual earith
*/
static inline void element_mul_mpz(element_t n, element_t a, mpz_t z) {
  PBC_ASSERT_MATCH2(n, a);
  n->field->mul_mpz(n, a, z);
}

/*@manual earith
Set 'n' = 'a' 'z', that is 'a' + 'a' + ... + 'a' where there are 'z' 'a'#'s#.
*/
static inline void element_mul_si(element_t n, element_t a,
                                  signed long int z) {
  PBC_ASSERT_MATCH2(n, a);
  n->field->mul_si(n, a, z);
}

/*@manual earith
'z' must be an element of a integer mod ring (i.e. *Z*~n~ for some n).
Set 'c' = 'a' 'z', that is 'a' + 'a' + ... + 'a'
where there are 'z' 'a''s.
*/
static inline void element_mul_zn(element_t c, element_t a, element_t z) {
  mpz_t z0;
  PBC_ASSERT_MATCH2(c, a);
  //TODO: check z->field is Zn
  mpz_init(z0);
  element_to_mpz(z0, z);
  element_mul_mpz(c, a, z0);
  mpz_clear(z0);
}

/*@manual earith
Set 'n' = 'a' / 'b'.
*/
static inline void element_div(element_t n, element_t a, element_t b) {
  PBC_ASSERT_MATCH3(n, a, b);
  n->field->div(n, a, b);
}

/*@manual earith
Set 'n' = 'a' + 'a'.
*/
static inline void element_double(element_t n, element_t a) {
  PBC_ASSERT_MATCH2(n, a);
  n->field->doub(n, a);
}

// Set n_i = a_i + a_i for all i at one time.
// Uses multi_doub(), which only elliptic curves have at the moment.
void element_multi_double(element_t n[], element_t a[], int m);

// Set n_i =a_i + b_i for all i at one time.
// Uses multi_add(), which only elliptic curves have at the moment.
void element_multi_add(element_t n[], element_t a[],element_t b[], int m);

/*@manual earith
Set 'n' = 'a/2'
*/
static inline void element_halve(element_t n, element_t a) {
  PBC_ASSERT_MATCH2(n, a);
  n->field->halve(n, a);
}

/*@manual earith
Set 'n' = 'a'^2^
*/
static inline void element_square(element_t n, element_t a) {
  PBC_ASSERT_MATCH2(n, a);
  n->field->square(n, a);
}

/*@manual epow
Set 'x' = 'a'^'n'^, that is
'a' times 'a' times ... times 'a' where there are 'n' 'a'#'s#.
*/
static inline void element_pow_mpz(element_t x, element_t a, mpz_t n) {
  PBC_ASSERT_MATCH2(x, a);
  x->field->pow_mpz(x, a, n);
}

/*@manual epow
Set 'x' = 'a'^'n'^, where 'n' is an element of a ring *Z*~N~
for some 'N' (typically the order of the algebraic structure 'x' lies in).
*/
static inline void element_pow_zn(element_t x, element_t a, element_t n) {
  mpz_t z;
  PBC_ASSERT_MATCH2(x, a);
  mpz_init(z);
  element_to_mpz(z, n);
  element_pow_mpz(x, a, z);
  mpz_clear(z);
}

/*@manual earith
Set 'n' = -'a'.
*/
static inline void element_neg(element_t n, element_t a) {
  PBC_ASSERT_MATCH2(n, a);
  n->field->neg(n, a);
}

/*@manual earith
Set 'n' to the inverse of 'a'.
*/
static inline void element_invert(element_t n, element_t a) {
  PBC_ASSERT_MATCH2(n, a);
  n->field->invert(n, a);
}

/*@manual erandom
If the 'e' lies in a finite algebraic structure,
assigns a uniformly random element to 'e'.
*/
static inline void element_random(element_t e) {
  e->field->random(e);
}

/*@manual ecmp
Returns true if 'n' is 1.
*/
static inline int element_is1(element_t n) {
  return n->field->is1(n);
}

/*@manual ecmp
Returns true if 'n' is 0.
*/
static inline int element_is0(element_t n) {
  return n->field->is0(n);
}

/*@manual ecmp
Returns 0 if 'a' and 'b' are the same, nonzero otherwise.
*/
static inline int element_cmp(element_t a, element_t b) {
  PBC_ASSERT_MATCH2(a, b);
  return a->field->cmp(a, b);
}

/*@manual ecmp
Returns nonzero if 'a' is a perfect square (quadratic residue),
zero otherwise.
*/
static inline int element_is_sqr(element_t a) {
  return a->field->is_sqr(a);
}

/*@manual ecmp
*/
static inline int element_sgn(element_t a) {
  return a->field->sign(a);
}

/*@manual ecmp
If 'a' is zero, returns 0. For nozero 'a' the behaviour depends on
the algebraic structure, but has the property that
element_sgn('a') = -element_sgn(-'a')
and
element_sgn('a') = 0 implies 'a' = 0 with overwhelming probability.
*/
static inline int element_sign(element_t a) {
  return a->field->sign(a);
}

static inline void element_sqrt(element_t a, element_t b) {
  PBC_ASSERT_MATCH2(a, b);
  a->field->sqrt(a, b);
}

/*@manual etrade
Returns the length in bytes the element 'e' will take to represent
*/
static inline int element_length_in_bytes(element_t e) {
  if (e->field->fixed_length_in_bytes < 0) {
    return e->field->length_in_bytes(e);
  } else {
    return e->field->fixed_length_in_bytes;
  }
}

/*@manual etrade
Converts 'e' to byte, writing the result in the buffer 'data'.
The number of bytes it will write can be determined from calling
*element_length_in_bytes()*. Returns number of bytes written.
*/
static inline int element_to_bytes(unsigned char *data, element_t e) {
  return e->field->to_bytes(data, e);
}

/*@manual etrade
Reads 'e' from the buffer 'data', and returns the number of bytes read.
*/
static inline int element_from_bytes(element_t e, unsigned char *data) {
  return e->field->from_bytes(e, data);
}

/*@manual epow
Sets 'x' = 'a1'^'n1'^ 'a2'^'n2'^, and is generally faster than
performing two separate exponentiations.
*/
void element_pow2_mpz(element_t x, element_t a1, mpz_t n1, element_t a2,
                      mpz_t n2);
/*@manual epow
Also sets 'x' = 'a1'^'n1'^ 'a2'^'n2'^,
but 'n1', 'n2' must be elements of a ring *Z*~n~ for some integer n.
*/
static inline void element_pow2_zn(element_t x, element_t a1, element_t n1,
                                   element_t a2, element_t n2) {
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
Sets 'x' = 'a1'^'n1'^ 'a2'^'n2'^ 'a3'^'n3'^,
generally faster than performing three separate exponentiations.
*/
void element_pow3_mpz(element_t x, element_t a1, mpz_t n1,
                      element_t a2, mpz_t n2, element_t a3, mpz_t n3);

/*@manual epow
Also sets 'x' = 'a1'^'n1'^ 'a2'^'n2'^ 'a3'^'n3'^,
but 'n1', 'n2', 'n3' must be elements of a ring *Z*~n~ for some integer n.
*/
static inline void element_pow3_zn(element_t x, element_t a1, element_t n1,
                                   element_t a2, element_t n2,
                                   element_t a3, element_t n3) {
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

static inline int mpz_is0(mpz_t z) {
  return !mpz_sgn(z);
  //return !mpz_cmp_ui(z, 0);
}

/*@manual etrade
Assumes 'e' is a point on an elliptic curve.
Writes the x-coordinate of 'e' to the buffer 'data'
*/
int element_to_bytes_x_only(unsigned char *data, element_t e);
/*@manual etrade
Assumes 'e' is a point on an elliptic curve.
Sets 'e' to a point with
x-coordinate represented by the buffer 'data'. This is not unique.
For each 'x'-coordinate, there exist two different points, at least
for the elliptic curves in PBC. (They are inverses of each other.)
*/
int element_from_bytes_x_only(element_t e, unsigned char *data);
/*@manual etrade
Assumes 'e' is a point on an elliptic curve.
Returns the length in bytes needed to hold the x-coordinate of 'e'.
*/
int element_length_in_bytes_x_only(element_t e);

/*@manual etrade
If possible, outputs a compressed form of the element 'e' to
the buffer of bytes 'data'.
Currently only implemented for points on an elliptic curve.
*/
int element_to_bytes_compressed(unsigned char *data, element_t e);

/*@manual etrade
Sets element 'e' to the element in compressed form in the buffer of bytes
'data'.
Currently only implemented for points on an elliptic curve.
*/
int element_from_bytes_compressed(element_t e, unsigned char *data);

/*@manual etrade
Returns the number of bytes needed to hold 'e' in compressed form.
Currently only implemented for points on an elliptic curve.
*/
int element_length_in_bytes_compressed(element_t e);

/*@manual epow
Prepare to exponentiate an element 'in', and store preprocessing information
in 'p'.
*/
static inline void element_pp_init(element_pp_t p, element_t in) {
  p->field = in->field;
  in->field->pp_init(p, in);
}

/*@manual epow
Clear 'p'. Should be called after 'p' is no longer needed.
*/
static inline void element_pp_clear(element_pp_t p) {
  p->field->pp_clear(p);
}

/*@manual epow
Raise 'in' to 'power' and store the result in 'out', where 'in'
is a previously preprocessed element, that is, the second argument
passed to a previous *element_pp_init* call.
*/
static inline void element_pp_pow(element_t out, mpz_ptr power,
                                  element_pp_t p) {
  p->field->pp_pow(out, power, p);
}

/*@manual epow
Same except 'power' is an element of *Z*~n~ for some integer n.
*/
static inline void element_pp_pow_zn(element_t out, element_t power,
                                     element_pp_t p) {
  mpz_t z;
  mpz_init(z);
  element_to_mpz(z, power);
  element_pp_pow(out, z, p);
  mpz_clear(z);
}

void pbc_mpz_out_raw_n(unsigned char *data, int n, mpz_t z);
void pbc_mpz_from_hash(mpz_t z, mpz_t limit,
                       unsigned char *data, unsigned int len);

/*@manual etrade
For points, returns the number of coordinates.
For polynomials, returns the number of coefficients.
Otherwise returns zero.
*/
static inline int element_item_count(element_t e) {
  return e->field->item_count(e);
}

/*@manual etrade
For points, returns 'n'#th# coordinate.
For polynomials, returns coefficient of 'x^n^'.
Otherwise returns NULL.
The element the return value points to may be modified.
*/
static inline element_ptr element_item(element_t e, int i) {
  // TODO: Document the following:
  // For polynomials, never zero the leading coefficient, e.g. never write:
  //  element_set0(element_item(f, poly_degree(f)));
  // Use poly_set_coeff0() to zero the leading coefficient.
  return e->field->item(e, i);
}

// Returns the field containing the items.
// Returns NULL if there are no items.
static inline field_ptr element_item_field(element_t e) {
  if (!element_item_count(e)) return NULL;
  return element_item(e, 0)->field;
}

/*@manual etrade
Equivalent to `element_item(a, 0)`.
*/
static inline element_ptr element_x(element_ptr a) {
  return a->field->get_x(a);
}
/*@manual etrade
Equivalent to `element_item(a, 1)`.
*/
static inline element_ptr element_y(element_ptr a) {
  return a->field->get_y(a);
}

/*@manual epow
Computes 'x' such that 'g^x^ = h' by brute force, where
'x' lies in a field where `element_set_mpz()` makes sense.
*/
void element_dlog_brute_force(element_t x, element_t g, element_t h);

/*@manual epow
Computes 'x' such that 'g^x^ = h' using Pollard rho method, where
'x' lies in a field where `element_set_mpz()` makes sense.
*/
void element_dlog_pollard_rho(element_t x, element_t g, element_t h);

// Trial division up to a given limit. If limit == NULL, then there is no limit.
// Call the callback for each factor found, abort and return 1 if the callback
// returns nonzero, otherwise return 0.
int pbc_trial_divide(int (*fun)(mpz_t factor,
                                unsigned int multiplicity,
                                void *scope_ptr),
                     void *scope_ptr,
                     mpz_t n,
                     mpz_ptr limit);

#endif  // __PBC_FIELD_H__
