/*
 * quadratic extension of any field
 */

//requires
// * field.h
#ifndef __PBC_FIELDQUADRATIC_H__
#define __PBC_FIELDQUADRATIC_H__

void field_init_quadratic(field_ptr f, field_ptr fbase);
void field_init_fi(field_ptr f, field_ptr fbase);

//naturally map an element from a field K to K[a]
void element_field_to_quadratic(element_ptr out, element_ptr in);

struct fq_data_s {
    //elements have the form x + ya
    //where a is the square root of some quadratic nonresidue in base field
    element_t x;
    element_t y;
};
typedef struct fq_data_s *fq_data_ptr;
typedef struct fq_data_s fq_data_t[1];

static inline element_ptr fi_re(element_ptr a)
{
    return ((fq_data_ptr) a->data)->x;
}

static inline element_ptr fi_im(element_ptr a)
{
    return ((fq_data_ptr) a->data)->y;
}

void element_field_to_fi(element_ptr a, element_ptr b);

#endif //__PBC_FIELDQUADRATIC_H__
