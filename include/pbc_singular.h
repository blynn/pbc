//requires
// * stdio.h
// * gmp.h
// * field.h
#ifndef SINGULAR_H
#define SINGULAR_H

void field_init_curve_singular_with_node(field_t c, field_t field);
void pairing_init_singular_with_node(pairing_t pairing, mpz_t q);

#endif //SINGULAR_H
