#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdint.h> // for intptr_t
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_field.h"
#include "pbc_multiz.h"
#include "pbc_poly.h"
#include "pbc_curve.h"
#include "pbc_memory.h"
#include "pbc_random.h"
#include "misc/darray.h"

// Per-field data.
typedef struct {
  field_ptr field; // The field where the curve is defined.
  element_t a, b;  // The curve is E: Y^2 = X^3 + a X + b.
  // cofac == NULL means we're using the whole group of points.
  // otherwise we're working in the subgroup of order #E / cofac,
  // where #E is the number of points in E.
  mpz_ptr cofac;
  // A generator of E.
  element_t gen_no_cofac;
  // A generator of the subgroup.
  element_t gen;
  // A non-NULL quotient_cmp means we are working with the quotient group of
  // order #E / quotient_cmp, and the points are actually coset
  // representatives. Thus for a comparison, we must multiply by quotient_cmp
  // before comparing.
  mpz_ptr quotient_cmp;
} *curve_data_ptr;

// Per-element data. Elements of this group are points on the elliptic curve.
typedef struct {
  int inf_flag;    // inf_flag == 1 means O, the point at infinity.
  element_t x, y;  // Otherwise we have the finite point (x, y).
} *point_ptr;

static void curve_init(element_ptr e) {
  curve_data_ptr cdp = e->field->data;
  point_ptr p = e->data = pbc_malloc(sizeof(*p));
  element_init(p->x, cdp->field);
  element_init(p->y, cdp->field);
  p->inf_flag = 1;
}

static void curve_clear(element_ptr e) {
  point_ptr p = e->data;
  element_clear(p->x);
  element_clear(p->y);
  pbc_free(e->data);
}

static int curve_is_valid_point(element_ptr e) {
  element_t t0, t1;
  int result;
  curve_data_ptr cdp = e->field->data;
  point_ptr p = e->data;

  if (p->inf_flag) return 1;

  element_init(t0, cdp->field);
  element_init(t1, cdp->field);
  element_square(t0, p->x);
  element_add(t0, t0, cdp->a);
  element_mul(t0, t0, p->x);
  element_add(t0, t0, cdp->b);
  element_square(t1, p->y);
  result = !element_cmp(t0, t1);

  element_clear(t0);
  element_clear(t1);
  return result;
}

static void curve_invert(element_ptr c, element_ptr a) {
  point_ptr r = c->data, p = a->data;

  if (p->inf_flag) {
    r->inf_flag = 1;
    return;
  }
  r->inf_flag = 0;
  element_set(r->x, p->x);
  element_neg(r->y, p->y);
}

static void curve_set(element_ptr c, element_ptr a) {
  point_ptr r = c->data, p = a->data;
  if (p->inf_flag) {
    r->inf_flag = 1;
    return;
  }
  r->inf_flag = 0;
  element_set(r->x, p->x);
  element_set(r->y, p->y);
}

static inline void double_no_check(point_ptr r, point_ptr p, element_ptr a) {
  element_t lambda, e0, e1;
  field_ptr f = r->x->field;

  element_init(lambda, f);
  element_init(e0, f);
  element_init(e1, f);

  //lambda = (3x^2 + a) / 2y
  element_square(lambda, p->x);
  element_mul_si(lambda, lambda, 3);
  element_add(lambda, lambda, a);

  element_double(e0, p->y);

  element_invert(e0, e0);
  element_mul(lambda, lambda, e0);
  //x1 = lambda^2 - 2x
  //element_add(e1, p->x, p->x);
  element_double(e1, p->x);
  element_square(e0, lambda);
  element_sub(e0, e0, e1);
  //y1 = (x - x1)lambda - y
  element_sub(e1, p->x, e0);
  element_mul(e1, e1, lambda);
  element_sub(e1, e1, p->y);

  element_set(r->x, e0);
  element_set(r->y, e1);
  r->inf_flag = 0;

  element_clear(lambda);
  element_clear(e0);
  element_clear(e1);
  return;
}

static void curve_double(element_ptr c, element_ptr a) {
  curve_data_ptr cdp = a->field->data;
  point_ptr r = c->data, p = a->data;
  if (p->inf_flag) {
    r->inf_flag = 1;
    return;
  }
  if (element_is0(p->y)) {
    r->inf_flag = 1;
    return;
  }
  double_no_check(r, p, cdp->a);
}

static void curve_mul(element_ptr c, element_ptr a, element_ptr b) {
  curve_data_ptr cdp = a->field->data;
  point_ptr r = c->data, p = a->data, q = b->data;

  if (p->inf_flag) {
    curve_set(c, b);
    return;
  }
  if (q->inf_flag) {
    curve_set(c, a);
    return;
  }
  if (!element_cmp(p->x, q->x)) {
    if (!element_cmp(p->y, q->y)) {
      if (element_is0(p->y)) {
        r->inf_flag = 1;
        return;
      } else {
        double_no_check(r, p, cdp->a);
        return;
      }
    }
    //points are inverses of each other
    r->inf_flag = 1;
    return;
  } else {
    element_t lambda, e0, e1;

    element_init(lambda, cdp->field);
    element_init(e0, cdp->field);
    element_init(e1, cdp->field);

    //lambda = (y2-y1)/(x2-x1)
    element_sub(e0, q->x, p->x);
    element_invert(e0, e0);
    element_sub(lambda, q->y, p->y);
    element_mul(lambda, lambda, e0);
    //x3 = lambda^2 - x1 - x2
    element_square(e0, lambda);
    element_sub(e0, e0, p->x);
    element_sub(e0, e0, q->x);
    //y3 = (x1-x3)lambda - y1
    element_sub(e1, p->x, e0);
    element_mul(e1, e1, lambda);
    element_sub(e1, e1, p->y);

    element_set(r->x, e0);
    element_set(r->y, e1);
    r->inf_flag = 0;

    element_clear(lambda);
    element_clear(e0);
    element_clear(e1);
  }
}

//compute c_i=a_i+a_i at one time.
static void multi_double(element_ptr c[], element_ptr a[], int n) {
  int i;
  element_t* table = pbc_malloc(sizeof(element_t)*n);  //a big problem?
  element_t e0, e1, e2;
  point_ptr q, r;
  curve_data_ptr cdp = a[0]->field->data;

  q=a[0]->data;
  element_init(e0,q->y->field);
  element_init(e1,q->y->field);
  element_init(e2,q->y->field);

  for(i=0; i<n; i++){
    q=a[i]->data; r=c[i]->data;
    element_init(table[i],q->y->field);

    if (q->inf_flag) {
      r->inf_flag = 1;
      continue;
    }
    if (element_is0(q->y)) {
      r->inf_flag = 1;
      continue;
    }
  }
  //to compute 1/2y multi. see Cohen's GTM139 Algorithm 10.3.4
  for(i=0; i<n; i++){
    q=a[i]->data;
    element_double(table[i],q->y);
    if(i>0) element_mul(table[i],table[i],table[i-1]);
  }
  element_invert(e2,table[n-1]); //ONLY ONE inv is required now.
  for(i=n-1; i>0; i--){
    q=a[i]->data;
    element_mul(table[i],table[i-1],e2);
    element_mul(e2,e2,q->y);
    element_double(e2,e2); //e2=e2*2y_j
  }
  element_set(table[0],e2); //e2 no longer used.

  for(i=0; i<n; i++){
    q=a[i]->data;
    r=c[i]->data;
    if(r->inf_flag) continue;

    //e2=lambda = (3x^2 + a) / 2y
    element_square(e2, q->x);
    element_mul_si(e2, e2, 3);
    element_add(e2, e2, cdp->a);

    element_mul(e2, e2, table[i]); //Recall that table[i]=1/2y_i
    //x1 = lambda^2 - 2x
    element_double(e1, q->x);
    element_square(e0, e2);
    element_sub(e0, e0, e1);
    //y1 = (x - x1)lambda - y
    element_sub(e1, q->x, e0);
    element_mul(e1, e1, e2);
    element_sub(e1, e1, q->y);
    element_set(r->x, e0);
    element_set(r->y, e1);
    r->inf_flag = 0;
  }

  element_clear(e0);
  element_clear(e1);
  element_clear(e2);
  for(i=0; i<n; i++){
    element_clear(table[i]);
  }
  pbc_free(table);
}

//compute c_i=a_i+b_i at one time.
static void multi_add(element_ptr c[], element_ptr a[], element_ptr b[], int n){
  int i;
  element_t* table = pbc_malloc(sizeof(element_t)*n);  //a big problem?
  point_ptr p, q, r;
  element_t e0, e1, e2;
  curve_data_ptr cdp = a[0]->field->data;

  p = a[0]->data;
  q = b[0]->data;
  element_init(e0, p->x->field);
  element_init(e1, p->x->field);
  element_init(e2, p->x->field);

  element_init(table[0], p->x->field);
  element_sub(table[0], q->x, p->x);
  for(i=1; i<n; i++){
    p = a[i]->data;
    q = b[i]->data;
    element_init(table[i], p->x->field);
    element_sub(table[i], q->x, p->x);
    element_mul(table[i], table[i], table[i-1]);
  }
  element_invert(e2, table[n-1]);
  for(i=n-1; i>0; i--){
    p = a[i]->data;
    q = b[i]->data;
    element_mul(table[i], table[i-1], e2);
    element_sub(e1, q->x, p->x);
    element_mul(e2,e2,e1); //e2=e2*(x2_j-x1_j)
  }
  element_set(table[0],e2); //e2 no longer used.

  for(i=0; i<n; i++){
    p = a[i]->data;
    q = b[i]->data;
    r = c[i]->data;
    if (p->inf_flag) {
      curve_set(c[i], b[i]);
      continue;
    }
    if (q->inf_flag) {
      curve_set(c[i], a[i]);
      continue;
    }
    if (!element_cmp(p->x, q->x)) { //a[i]=b[i]
      if (!element_cmp(p->y, q->y)) {
        if (element_is0(p->y)) {
          r->inf_flag = 1;
          continue;
        } else {
          double_no_check(r, p, cdp->a);
          continue;
        }
      }
      //points are inverses of each other
      r->inf_flag = 1;
      continue;
    } else {
      //lambda = (y2-y1)/(x2-x1)
      element_sub(e2, q->y, p->y);
      element_mul(e2, e2, table[i]);
      //x3 = lambda^2 - x1 - x2
      element_square(e0, e2);
      element_sub(e0, e0, p->x);
      element_sub(e0, e0, q->x);
      //y3 = (x1-x3)lambda - y1
      element_sub(e1, p->x, e0);
      element_mul(e1, e1, e2);
      element_sub(e1, e1, p->y);
      element_set(r->x, e0);
      element_set(r->y, e1);
      r->inf_flag = 0;
    }
  }
  element_clear(e0);
  element_clear(e1);
  element_clear(e2);
  for(i=0; i<n; i++){
    element_clear(table[i]);
  }
  pbc_free(table);
}


static inline int point_cmp(point_ptr p, point_ptr q) {
  if (p->inf_flag || q->inf_flag) {
    return !(p->inf_flag && q->inf_flag);
  }
  return element_cmp(p->x, q->x) || element_cmp(p->y, q->y);
}

static int curve_cmp(element_ptr a, element_ptr b) {
  if (a == b) {
    return 0;
  } else {
    // If we're working with a quotient group we must account for different
    // representatives of the same coset.
    curve_data_ptr cdp = a->field->data;
    if (cdp->quotient_cmp) {
      element_t e;
      element_init_same_as(e, a);
      element_div(e, a, b);
      element_pow_mpz(e, e, cdp->quotient_cmp);
      int result = !element_is1(e);
      element_clear(e);
      return result;
    }
    return point_cmp(a->data, b->data);
  }
}

static void curve_set1(element_ptr x) {
  point_ptr p = x->data;
  p->inf_flag = 1;
}

static int curve_is1(element_ptr x) {
  point_ptr p = x->data;
  return p->inf_flag;
}

static void curve_random_no_cofac_solvefory(element_ptr a) {
  //TODO: with 0.5 probability negate y-coord
  curve_data_ptr cdp = a->field->data;
  point_ptr p = a->data;
  element_t t;

  element_init(t, cdp->field);
  p->inf_flag = 0;
  do {
    element_random(p->x);
    element_square(t, p->x);
    element_add(t, t, cdp->a);
    element_mul(t, t, p->x);
    element_add(t, t, cdp->b);
  } while (!element_is_sqr(t));
  element_sqrt(p->y, t);
  element_clear(t);
}

static void curve_random_solvefory(element_ptr a) {
  curve_data_ptr cdp = a->field->data;
  curve_random_no_cofac_solvefory(a);
  if (cdp->cofac) element_mul_mpz(a, a, cdp->cofac);
}

static void curve_random_pointmul(element_ptr a) {
  curve_data_ptr cdp = a->field->data;
  mpz_t x;
  mpz_init(x);

  pbc_mpz_random(x, a->field->order);
  element_mul_mpz(a, cdp->gen, x);
  mpz_clear(x);
}

void field_curve_use_random_solvefory(field_ptr f) {
  f->random = curve_random_solvefory;
}

void curve_set_gen_no_cofac(element_ptr a) {
  curve_data_ptr cdp = a->field->data;
  element_set(a, cdp->gen_no_cofac);
}

static int curve_sign(element_ptr e) {
  point_ptr p = e->data;
  if (p->inf_flag) return 0;
  return element_sign(p->y);
}

static void curve_from_hash(element_t a, void *data, int len) {
  element_t t, t1;
  point_ptr p = a->data;
  curve_data_ptr cdp = a->field->data;

  element_init(t, cdp->field);
  element_init(t1, cdp->field);
  p->inf_flag = 0;
  element_from_hash(p->x, data, len);
  for(;;) {
    element_square(t, p->x);
    element_add(t, t, cdp->a);
    element_mul(t, t, p->x);
    element_add(t, t, cdp->b);
    if (element_is_sqr(t)) break;
    // Compute x <- x^2 + 1 and try again.
    element_square(p->x, p->x);
    element_set1(t);
    element_add(p->x, p->x, t);
  }
  element_sqrt(p->y, t);
  if (element_sgn(p->y) < 0) element_neg(p->y, p->y);

  if (cdp->cofac) element_mul_mpz(a, a, cdp->cofac);

  element_clear(t);
  element_clear(t1);
}

static size_t curve_out_str(FILE *stream, int base, element_ptr a) {
  point_ptr p = a->data;
  size_t result, status;
  if (p->inf_flag) {
    if (EOF == fputc('O', stream)) return 0;
    return 1;
  }
  if (EOF == fputc('[', stream)) return 0;
  result = element_out_str(stream, base, p->x);
  if (!result) return 0;
  if (EOF == fputs(", ", stream)) return 0;
  status = element_out_str(stream, base, p->y);
  if (!status) return 0;
  if (EOF == fputc(']', stream)) return 0;
  return result + status + 4;
}

static int curve_snprint(char *s, size_t n, element_ptr a) {
  point_ptr p = a->data;
  size_t result = 0, left;
  int status;

  #define clip_sub() {                   \
    result += status;                    \
    left = result >= n ? 0 : n - result; \
  }

  if (p->inf_flag) {
    status = snprintf(s, n, "O");
    if (status < 0) return status;
    return 1;
  }

  status = snprintf(s, n, "[");
  if (status < 0) return status;
  clip_sub();
  status = element_snprint(s + result, left, p->x);
  if (status < 0) return status;
  clip_sub();
  status = snprintf(s + result, left, ", ");
  if (status < 0) return status;
  clip_sub();
  status = element_snprint(s + result, left, p->y);
  if (status < 0) return status;
  clip_sub();
  status = snprintf(s + result, left, "]");
  if (status < 0) return status;
  return result + status;
  #undef clip_sub
}

static void curve_set_multiz(element_ptr a, multiz m) {
  if (multiz_is_z(m)) {
    if (multiz_is0(m)) {
      element_set0(a);
      return;
    }
    pbc_warn("bad multiz");
    return;
  } else {
    if (multiz_count(m) < 2) {
      pbc_warn("multiz has too few coefficients");
      return;
    }
    point_ptr p = a->data;
    p->inf_flag = 0;
    element_set_multiz(p->x, multiz_at(m, 0));
    element_set_multiz(p->y, multiz_at(m, 1));
  }
}

static int curve_set_str(element_ptr e, const char *s, int base) {
  point_ptr p = e->data;
  const char *cp = s;
  element_set0(e);
  while (*cp && isspace(*cp)) cp++;
  if (*cp == 'O') {
    return cp - s + 1;
  }
  p->inf_flag = 0;
  if (*cp != '[') return 0;
  cp++;
  cp += element_set_str(p->x, cp, base);
  while (*cp && isspace(*cp)) cp++;
  if (*cp != ',') return 0;
  cp++;
  cp += element_set_str(p->y, cp, base);
  if (*cp != ']') return 0;

  if (!curve_is_valid_point(e)) {
    element_set0(e);
    return 0;
  }
  return cp - s + 1;
}

static void field_clear_curve(field_t f) {
  curve_data_ptr cdp;
  cdp = f->data;
  element_clear(cdp->gen);
  element_clear(cdp->gen_no_cofac);
  if (cdp->cofac) {
    mpz_clear(cdp->cofac);
    pbc_free(cdp->cofac);
  }
  if (cdp->quotient_cmp) {
    mpz_clear(cdp->quotient_cmp);
    pbc_free(cdp->quotient_cmp);
  }
  element_clear(cdp->a);
  element_clear(cdp->b);
  pbc_free(cdp);
}

static int curve_length_in_bytes(element_ptr x) {
  point_ptr p = x->data;
  return element_length_in_bytes(p->x) + element_length_in_bytes(p->y);
}

static int curve_to_bytes(unsigned char *data, element_t e) {
  point_ptr P = e->data;
  int len;
  len = element_to_bytes(data, P->x);
  len += element_to_bytes(data + len, P->y);
  return len;
}

static int curve_from_bytes(element_t e, unsigned char *data) {
  point_ptr P = e->data;
  int len;

  P->inf_flag = 0;
  len = element_from_bytes(P->x, data);
  len += element_from_bytes(P->y, data + len);
  //if point does not lie on curve, set it to O
  if (!curve_is_valid_point(e)) {
    element_set0(e);
  }
  return len;
}

static void curve_out_info(FILE *out, field_t f) {
  int len;
  fprintf(out, "elliptic curve");
  if ((len = f->fixed_length_in_bytes)) {
    fprintf(out, ", bits per coord = %d", len * 8 / 2);
  } else {
    fprintf(out, "variable-length");
  }
}

static int odd_curve_is_sqr(element_ptr e) {
  UNUSED_VAR(e);
  return 1;
}

//TODO: untested
static int even_curve_is_sqr(element_ptr e) {
  mpz_t z;
  element_t e1;
  int result;

  mpz_init(z);
  element_init(e1, e->field);
  mpz_sub_ui(z, e->field->order, 1);
  mpz_fdiv_q_2exp(z, z, 1);
  element_pow_mpz(e1, e, z);
  result = element_is1(e1);

  mpz_clear(z);
  element_clear(e1);
  return result;
}

static int curve_item_count(element_ptr e) {
  if (element_is0(e)) {
    return 0;
  }
  return 2;
}

static element_ptr curve_item(element_ptr e, int i) {
  if (element_is0(e)) return NULL;
  point_ptr P = e->data;
  switch(i) {
    case 0:
      return P->x;
    case 1:
      return P->y;
    default:
      return NULL;
  }
}

static element_ptr curve_get_x(element_ptr e) {
  point_ptr P = e->data;
  return P->x;
}

static element_ptr curve_get_y(element_ptr e) {
  point_ptr P = e->data;
  return P->y;
}

void field_init_curve_ab(field_ptr f, element_ptr a, element_ptr b, mpz_t order, mpz_t cofac) {
  /*
  if (element_is0(a)) {
    c->double_nocheck = cc_double_no_check_ais0;
  } else {
    c->double_nocheck = cc_double_no_check;
  }
  */
  curve_data_ptr cdp;
  field_init(f);
  mpz_set(f->order, order);
  cdp = f->data = pbc_malloc(sizeof(*cdp));
  cdp->field = a->field;
  element_init(cdp->a, cdp->field);
  element_init(cdp->b, cdp->field);
  element_set(cdp->a, a);
  element_set(cdp->b, b);

  f->init = curve_init;
  f->clear = curve_clear;
  f->neg = f->invert = curve_invert;
  f->square = f->doub = curve_double;
  f->multi_doub = multi_double;
  f->add = f->mul = curve_mul;
  f->multi_add = multi_add;
  f->mul_mpz = element_pow_mpz;
  f->cmp = curve_cmp;
  f->set0 = f->set1 = curve_set1;
  f->is0 = f->is1 = curve_is1;
  f->sign = curve_sign;
  f->set = curve_set;
  f->random = curve_random_pointmul;
  //f->random = curve_random_solvefory;
  f->from_hash = curve_from_hash;
  f->out_str = curve_out_str;
  f->snprint = curve_snprint;
  f->set_multiz = curve_set_multiz;
  f->set_str = curve_set_str;
  f->field_clear = field_clear_curve;
  if (cdp->field->fixed_length_in_bytes < 0) {
    f->length_in_bytes = curve_length_in_bytes;
  } else {
    f->fixed_length_in_bytes = 2 * cdp->field->fixed_length_in_bytes;
  }
  f->to_bytes = curve_to_bytes;
  f->from_bytes = curve_from_bytes;
  f->out_info = curve_out_info;
  f->item_count = curve_item_count;
  f->item = curve_item;
  f->get_x = curve_get_x;
  f->get_y = curve_get_y;

  if (mpz_odd_p(order)) {
    f->is_sqr = odd_curve_is_sqr;
  } else {
    f->is_sqr = even_curve_is_sqr;
  }

  element_init(cdp->gen_no_cofac, f);
  element_init(cdp->gen, f);
  curve_random_no_cofac_solvefory(cdp->gen_no_cofac);
  if (cofac) {
    cdp->cofac = pbc_malloc(sizeof(mpz_t));
    mpz_init(cdp->cofac);
    mpz_set(cdp->cofac, cofac);
    element_mul_mpz(cdp->gen, cdp->gen_no_cofac, cofac);
  } else{
    cdp->cofac = NULL;
    element_set(cdp->gen, cdp->gen_no_cofac);
  }
  cdp->quotient_cmp = NULL;
}

// Requires e to be a point on an elliptic curve.
int element_to_bytes_compressed(unsigned char *data, element_ptr e) {
  point_ptr P = e->data;
  int len;
  len = element_to_bytes(data, P->x);
  if (element_sign(P->y) > 0) {
    data[len] = 1;
  } else {
    data[len] = 0;
  }
  len++;
  return len;
}

// Computes a point on the elliptic curve Y^2 = X^3 + a X + b given its
// x-coordinate.
// Requires a solution to exist.
static void point_from_x(point_ptr p, element_t x, element_t a, element_t b) {
  element_t t;

  element_init(t, x->field);
  p->inf_flag = 0;
  element_square(t, x);
  element_add(t, t, a);
  element_mul(t, t, x);
  element_add(t, t, b);
  element_sqrt(p->y, t);
  element_set(p->x, x);

  element_clear(t);
}

void curve_from_x(element_ptr e, element_t x) {
  curve_data_ptr cdp = e->field->data;
  point_from_x(e->data, x, cdp->a, cdp->b);
}

// Requires e to be a point on an elliptic curve.
int element_from_bytes_compressed(element_ptr e, unsigned char *data) {
  curve_data_ptr cdp = e->field->data;
  point_ptr P = e->data;
  int len;
  len = element_from_bytes(P->x, data);
  point_from_x(P, P->x, cdp->a, cdp->b);

  if (data[len]) {
    if (element_sign(P->y) < 0) element_neg(P->y, P->y);
  } else if (element_sign(P->y) > 0) {
    element_neg(P->y, P->y);
  }
  len++;
  return len;
}

int element_length_in_bytes_compressed(element_ptr e) {
  point_ptr P = e->data;
  return element_length_in_bytes(P->x) + 1;
}

// Requires e to be a point on an elliptic curve.
int element_to_bytes_x_only(unsigned char *data, element_ptr e) {
  point_ptr P = e->data;
  int len;
  len = element_to_bytes(data, P->x);
  return len;
}

// Requires e to be a point on an elliptic curve.
int element_from_bytes_x_only(element_ptr e, unsigned char *data) {
  curve_data_ptr cdp = e->field->data;
  point_ptr P = e->data;
  int len;
  len = element_from_bytes(P->x, data);
  point_from_x(P, P->x, cdp->a, cdp->b);
  return len;
}

int element_length_in_bytes_x_only(element_ptr e) {
  point_ptr P = e->data;
  return element_length_in_bytes(P->x);
}

inline element_ptr curve_x_coord(element_t e) {
  return ((point_ptr) e->data)->x;
}

inline element_ptr curve_y_coord(element_t e) {
  return ((point_ptr) e->data)->y;
}

inline element_ptr curve_a_coeff(element_t e) {
  return ((curve_data_ptr) e->field->data)->a;
}

inline element_ptr curve_b_coeff(element_t e) {
  return ((curve_data_ptr) e->field->data)->b;
}

inline element_ptr curve_field_a_coeff(field_t f) {
  return ((curve_data_ptr) f->data)->a;
}

inline element_ptr curve_field_b_coeff(field_t f) {
  return ((curve_data_ptr) f->data)->b;
}

void field_init_curve_ab_map(field_t cnew, field_t c,
    fieldmap map, field_ptr mapdest,
    mpz_t ordernew, mpz_t cofacnew) {
  element_t a, b;
  curve_data_ptr cdp = c->data;

  element_init(a, mapdest);
  element_init(b, mapdest);

  map(a, cdp->a);
  map(b, cdp->b);

  field_init_curve_ab(cnew, a, b, ordernew, cofacnew);
  element_clear(a);
  element_clear(b);
}

// Existing points are invalidated as this mangles c.
void field_reinit_curve_twist(field_ptr c) {
  curve_data_ptr cdp = c->data;
  element_ptr nqr = field_get_nqr(cdp->field);
  element_mul(cdp->a, cdp->a, nqr);
  element_mul(cdp->a, cdp->a, nqr);
  element_mul(cdp->b, cdp->b, nqr);
  element_mul(cdp->b, cdp->b, nqr);
  element_mul(cdp->b, cdp->b, nqr);

  // Recompute generators.
  curve_random_no_cofac_solvefory(cdp->gen_no_cofac);
  if (cdp->cofac) {
    element_mul_mpz(cdp->gen, cdp->gen_no_cofac, cdp->cofac);
  } else{
    element_set(cdp->gen, cdp->gen_no_cofac);
  }
}

// I could generalize this for all fields, but is there any point?
void field_curve_set_quotient_cmp(field_ptr c, mpz_t quotient_cmp) {
  curve_data_ptr cdp = c->data;
  cdp->quotient_cmp = pbc_malloc(sizeof(mpz_t));
  mpz_init(cdp->quotient_cmp);
  mpz_set(cdp->quotient_cmp, quotient_cmp);
}

// Requires j != 0, 1728.
void field_init_curve_j(field_ptr f, element_ptr j, mpz_t order, mpz_t cofac) {
  element_t a, b;
  element_init(a, j->field);
  element_init(b, j->field);

  element_set_si(a, 1728);
  element_sub(a, a, j);
  element_invert(a, a);
  element_mul(a, a, j);

  //b = 2 j / (1728 - j)
  element_add(b, a, a);
  //a = 3 j / (1728 - j)
  element_add(a, a, b);
  field_init_curve_ab(f, a, b, order, cofac);

  element_clear(a);
  element_clear(b);
}

void field_init_curve_b(field_ptr f, element_ptr b, mpz_t order, mpz_t cofac) {
  element_t a;
  element_init(a, b->field);
  field_init_curve_ab(f, a, b, order, cofac);

  element_clear(a);
}

// Compute trace of Frobenius at q^n given trace at q.
// See p.105 of Blake, Seroussi and Smart.
void pbc_mpz_trace_n(mpz_t res, mpz_t q, mpz_t trace, int n) {
  int i;
  mpz_t c0, c1, c2;
  mpz_t t0;

  mpz_init(c0);
  mpz_init(c1);
  mpz_init(c2);
  mpz_init(t0);
  mpz_set_ui(c2, 2);
  mpz_set(c1, trace);
  for (i=2; i<=n; i++) {
    mpz_mul(c0, trace, c1);
    mpz_mul(t0, q, c2);
    mpz_sub(c0, c0, t0);
    mpz_set(c2, c1);
    mpz_set(c1, c0);
  }
  mpz_set(res, c1);
  mpz_clear(t0);
  mpz_clear(c2);
  mpz_clear(c1);
  mpz_clear(c0);
}

// Given q, t such that #E(F_q) = q - t + 1, compute #E(F_q^k).
void pbc_mpz_curve_order_extn(mpz_t res, mpz_t q, mpz_t t, int k) {
  mpz_t z;
  mpz_t tk;
  mpz_init(z);
  mpz_init(tk);
  mpz_pow_ui(z, q, k);
  mpz_add_ui(z, z, 1);
  pbc_mpz_trace_n(tk, q, t, k);
  mpz_sub(z, z, tk);
  mpz_set(res, z);
  mpz_clear(z);
  mpz_clear(tk);
}

void curve_set_si(element_t R, long int x, long int y) {
  point_ptr p = R->data;
  element_set_si(p->x, x);
  element_set_si(p->y, y);
  p->inf_flag = 0;
}
