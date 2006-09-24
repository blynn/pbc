//requires
// * stdio.h
// * gmp.h
// * fops.h
// * pairing.h
#ifndef F_PARAM_H
#define F_PARAM_H

struct f_param_s {
    mpz_t q; //curve defined over F_q
    mpz_t r; //r is the order of the curve
    mpz_t b; //curve equation is y^2 = x^3 + b
    mpz_t beta; //beta is a quadratic nonresidue in Fq
	//we use F_q^2 = F_q[sqrt(beta)]
    mpz_t alpha0, alpha1;
	//the polynomial x^6 + alpha0 + alpha1 sqrt(beta)
	//is irreducible over F_q^2[x], so
	//we can extend F_q^2 to F_q^12 using the
	//sixth root of -(alpha0 + alpha1 sqrt(beta))
};
typedef struct f_param_s f_param_t[1];
typedef struct f_param_s *f_param_ptr;

/*@manual fparam
Initialize ''p''. This must be called before ''p'' can be used.
*/
void f_param_init(f_param_t fp);

/*@manual fparam
Clear ''p''. This should be called after ''p'' is no longer needed.
*/
void f_param_clear(f_param_t fp);

/*@manual fparam
Generate type F pairing parameters and store them in ''p''.
Both the group order r and the order of the base field q will be roughly
''bits''-bit numbers.
To be secure, generic discrete log algorithms must
be infeasible in groups of order r, and finite field discrete log algorithms
must be infeasible in finite fields of order q^12.
Typical value: ''bits'' = 160.
*/
void f_param_gen(f_param_t fp, int bits);

/*@manual fparam
Write the parameters in ''p'' in a text format onto ''stream''.
*/
void f_param_out_str(FILE *stream, f_param_ptr p);
void f_param_inp_generic (f_param_ptr p, fetch_ops_t fops, void *ctx);
void pairing_init_f_param(pairing_t pairing, f_param_t param);

#endif //F_PARAM_H
