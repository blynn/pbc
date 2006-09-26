/* There does not appear to be a succint name for rings of type Z/nZ.
 * Sage calls it integer mod ring.
 * NTL calls it ZZ_p.
 * I'll call it fp, as it's the quickest to type.
 * "zn" might be better since it can also handle composite numbers.
 */
// Requires
// * field.h
#ifndef FP_H
#define FP_H

void field_init_naive_fp(field_ptr f, mpz_t prime);
void field_init_fast_fp(field_ptr f, mpz_t prime);

static inline void field_init_fp(field_ptr f, mpz_t prime)
{
    //field_init_naive_fp(f, prime);
    field_init_fast_fp(f, prime);
}

#endif //FP_H
