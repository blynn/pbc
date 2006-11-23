// requires
// * poly.h
#ifndef __PBC_HILBERT_H__
#define __PBC_HILBERT_H__

void hilbert_poly(darray_t P, int D);
void hilbert_poly_clear(darray_t P);
int findroot(element_ptr root, element_ptr poly);

#endif //__PBC_HILBERT_H__
