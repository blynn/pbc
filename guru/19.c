/*
 * Toy example of a field where the Tate pairing can be used
 * but the Weil pairing cannot.
 *
 * Consider the curve E: y^2 = x^3 + x + 6 over F_19:
 * E(F_19) is a cyclic group of order 18.
 * Thus E[3] is not contained in F_19
 * (it turns out E[3] is contained in F_19^3).
 *
 * Hence the Weil pairing cannot be defined over F_19
 * However, F_19 contains the cube roots of unity
 * so we can compute the Tate pairing
 */

/*
 * P = (12,13) generates a group of order 3:
 * <(12,13)> = {(12,13), (12,6), O}
 * e(P,P) = 7, so we have the isomorphism
 * <(12,13)> = <7> (in F_19^*)
 *
 * Similarly P = (4, 6) generates a group of order 9, and we find
 * <(4,6)> = <4>
 *
 * P = (0, 5) generates all of E(F_19)
 * Miller's algorithm will not allow us to calculate e(P, P) without
 * first extending F_19.
 * Instead of extending, we could manipulate rational functions since
 * 19 is small enough that an explicit expression of f_P can be found.
 */

#include "pbc.h"
#include "pbc_fp.h"
#include "pbc_fieldquadratic.h"

static void miller(element_t res, element_t P, element_ptr QR, element_ptr R, int n) {
  // Collate divisions.
  int m;
  element_t v, vd;
  element_t Z;
  element_t a, b, c;
  const element_ptr cca = curve_a_coeff(P);
  const element_ptr Px = curve_x_coord(P);
  const element_ptr Py = curve_y_coord(P);
  element_t e0, e1;
  mpz_t q;
  element_ptr Zx, Zy;
  const element_ptr numx = curve_x_coord(QR);
  const element_ptr numy = curve_y_coord(QR);
  const element_ptr denomx = curve_x_coord(R);
  const element_ptr denomy = curve_y_coord(R);

  void do_vertical(element_t e, element_t edenom)
  {
    element_sub(e0, numx, Zx);
    element_mul(e, e, e0);

    element_sub(e0, denomx, Zx);
    element_mul(edenom, edenom, e0);
  }

  void do_tangent(element_t e, element_t edenom)
  {
    //a = -slope_tangent(A.x, A.y);
    //b = 1;
    //c = -(A.y + a * A.x);
    //but we multiply by 2*A.y to avoid division

    //a = -Ax * (Ax + Ax + Ax + twicea_2) - a_4;
    //Common curves: a2 = 0 (and cc->a is a_4), so
    //a = -(3 Ax^2 + cc->a)
    //b = 2 * Ay
    //c = -(2 Ay^2 + a Ax);

    if (element_is0(Zy)) {
      do_vertical(e, edenom);
      return;
    }
    element_square(a, Zx);
    element_mul_si(a, a, 3);
    element_add(a, a, cca);
    element_neg(a, a);

    element_add(b, Zy, Zy);

    element_mul(e0, b, Zy);
    element_mul(c, a, Zx);
    element_add(c, c, e0);
    element_neg(c, c);

    element_mul(e0, a, numx);
    element_mul(e1, b, numy);
    element_add(e0, e0, e1);
    element_add(e0, e0, c);
    element_mul(e, e, e0);

    element_mul(e0, a, denomx);
    element_mul(e1, b, denomy);
    element_add(e0, e0, e1);
    element_add(e0, e0, c);
    element_mul(edenom, edenom, e0);
  }

  void do_line(element_ptr e, element_ptr edenom)
  {
    if (!element_cmp(Zx, Px)) {
      if (!element_cmp(Zy, Py)) {
        do_tangent(e, edenom);
      } else {
        do_vertical(e, edenom);
      }
      return;
    }

    element_sub(b, Px, Zx);
    element_sub(a, Zy, Py);
    element_mul(c, Zx, Py);
    element_mul(e0, Zy, Px);
    element_sub(c, c, e0);

    element_mul(e0, a, numx);
    element_mul(e1, b, numy);
    element_add(e0, e0, e1);
    element_add(e0, e0, c);
    element_mul(e, e, e0);

    element_mul(e0, a, denomx);
    element_mul(e1, b, denomy);
    element_add(e0, e0, e1);
    element_add(e0, e0, c);
    element_mul(edenom, edenom, e0);
  }

  element_init(a, res->field);
  element_init(b, res->field);
  element_init(c, res->field);
  element_init(e0, res->field);
  element_init(e1, res->field);

  element_init(v, res->field);
  element_init(vd, res->field);
  element_init(Z, P->field);

  element_set(Z, P);
  Zx = curve_x_coord(Z);
  Zy = curve_y_coord(Z);

  element_set1(v);
  element_set1(vd);

  mpz_init(q);
  mpz_set_ui(q, n);
  m = mpz_sizeinbase(q, 2) - 2;

  while(m >= 0) {
    element_square(v, v);
    element_square(vd, vd);
    do_tangent(v, vd);
    element_double(Z, Z);
    do_vertical(vd, v);

    if (mpz_tstbit(q, m)) {
      do_line(v, vd);
      element_add(Z, Z, P);
      if (m) {
        do_vertical(vd, v);
      }
    }
    m--;
  }

  mpz_clear(q);

  element_invert(vd, vd);
  element_mul(res, v, vd);

  element_clear(v);
  element_clear(vd);
  element_clear(Z);
  element_clear(a);
  element_clear(b);
  element_clear(c);
  element_clear(e0);
  element_clear(e1);
}

static void tate_3(element_ptr out, element_ptr P, element_ptr Q, element_ptr R) {
  mpz_t six;

  mpz_init(six);
  mpz_set_ui(six, 6);
  element_t QR;
  element_t e0;

  element_init(QR, P->field);
  element_init(e0, out->field);

  element_add(QR, Q, R);

  //for subgroup size 3, -2P = P, hence
  //the tangent line at P has divisor 3(P) - 3(O)

  miller(out, P, QR, R, 3);

  element_pow_mpz(out, out, six);
  element_clear(QR);
  element_clear(e0);
  mpz_clear(six);
}

static void tate_9(element_ptr out, element_ptr P, element_ptr Q, element_ptr R) {
  element_t QR;
  element_init(QR, P->field);

  element_add(QR, Q, R);

  miller(out, P, QR, R, 9);

  element_square(out, out);

  element_clear(QR);
}

static void tate_18(element_ptr out, element_ptr P, element_ptr Q, element_ptr R, element_ptr S) {
  mpz_t pow;
  element_t PR;
  element_t QS;
  element_init(PR, P->field);
  element_init(QS, P->field);
  element_t outd;

  element_init(outd, out->field);

  mpz_init(pow);
  mpz_set_ui(pow, (19*19-1)/18);

  element_add(PR, P, R);
  element_add(QS, Q, S);

  if (element_is0(QS)) {
    element_t S2;
    element_init(S2, P->field);
    element_double(S2, S);
    miller(out, PR, S, S2, 18);
    miller(outd, R, S, S2, 18);
    element_clear(S2);
  } else {
    miller(out, PR, QS, S, 18);
    miller(outd, R, QS, S, 18);
  }

  element_clear(PR);
  element_clear(QS);

  element_invert(outd, outd);
  element_mul(out, out, outd);
  element_pow_mpz(out, out, pow);

  element_clear(outd);
  mpz_clear(pow);
}

int main(void) {
  field_t c;
  field_t Z19;
  element_t P, Q, R;
  mpz_t q, z;
  element_t a, b;
  int i;

  field_t Z19_2;
  field_t c2;
  element_t P2, Q2, R2;
  element_t a2;

  mpz_init(q);
  mpz_init(z);

  mpz_set_ui(q, 19);

  field_init_fp(Z19, q);
  element_init(a, Z19);
  element_init(b, Z19);

  element_set_si(a, 1);
  element_set_si(b, 6);

  mpz_set_ui(q, 18);
  field_init_curve_ab(c, a, b, q, NULL);
  element_init(P, c);
  element_init(Q, c);
  element_init(R, c);

  printf("Y^2 = X^3 + X + 6 over F_19\n");
  //(0,+/-5) is a generator
  element_set0(a);
  curve_from_x(R, a);

  for (i=1; i<19; i++) {
    mpz_set_si(z, i);
    element_mul_mpz(Q, R, z);
    element_printf("%dR = %B\n", i, Q);
  }

  mpz_set_ui(z, 6);
  element_mul_mpz(P, R, z);
  //P has order 3
  element_printf("P = %B\n", P);

  for (i=1; i<=3; i++) {
    mpz_set_si(z, i);
    element_mul_mpz(Q, R, z);
    tate_3(a, P, Q, R);
    element_printf("e_3(P,%dR) = %B\n", i, a);
  }

  element_double(P, R);
  //P has order 9
  element_printf("P = %B\n", P);
  for (i=1; i<=9; i++) {
    mpz_set_si(z, i);
    //we're supposed to use multiples of R
    //but 2R works just as well and it allows us
    //to use R as the offset every time
    element_mul_mpz(Q, P, z);
    tate_9(a, P, Q, R);
    element_printf("e_9(P,%dP) = %B\n", i, a);
  }

  //to do the pairing on all of E(F_19) we need to move to F_19^2
  //or compute the rational function explicitly
  printf("moving to F_19^2\n");
  field_init_fi(Z19_2, Z19);

  //don't need to tell it the real order
  field_init_curve_ab_map(c2, c, element_field_to_fi, Z19_2, q, NULL);
  element_init(P2, c2);
  element_init(Q2, c2);
  element_init(R2, c2);

  element_init(a2, Z19_2);
  element_set0(a2);
  curve_from_x(P2, a2);

  element_random(R2);

  element_printf("P = %B\n", P2);

  for (i=1; i<=18; i++) {
    mpz_set_si(z, i);
    element_mul_mpz(Q2, P2, z);
    tate_18(a2, P2, Q2, R2, P2);
    element_printf("e_18(P,%dP) = %B\n", i, a2);
  }

  element_clear(P2);
  element_clear(Q2);
  element_clear(R2);
  element_clear(a2);
  field_clear(c2);
  field_clear(Z19_2);

  field_clear(c);
  element_clear(a);
  element_clear(b);
  element_clear(P);
  element_clear(Q);
  element_clear(R);
  field_clear(Z19);

  mpz_clear(q);
  mpz_clear(z);
  return 0;
}
