// requires
// * poly.h
#ifndef HILBERT_H
#define HILBERT_H

void hilbert_poly(darray_t P, int D);
void hilbert_poly_clear(darray_t P);
int findroot(element_ptr root, element_ptr poly);

#endif //HILBERT_H
