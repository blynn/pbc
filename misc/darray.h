// darray = "dynamic array"
// A linked-list implementation using C arrays.

#ifndef __PBC_DARRAY_H__
#define __PBC_DARRAY_H__

#pragma GCC visibility push(hidden)

struct darray_s {
  void **item;
  int count;
  int max;
};

typedef struct darray_s darray_t[1];
typedef struct darray_s *darray_ptr;

/*@manual darray
Initialize a dynamic array 'a'. Must be called before 'a' is used.
*/
void darray_init(darray_t a);
darray_ptr darray_new(void);

void darray_free(darray_ptr a);

/*@manual darray
Clears a dynamic array 'a'. Should be called after 'a' is no longer needed.
*/
void darray_clear(darray_t a);

/*@manual darray
Appends 'p' to the dynamic array 'a'.
*/
void darray_append(darray_t a, void *p);

/*@manual darray
Returns the pointer at index 'i' in the dynamic array 'a'.
*/
static inline void *darray_at(darray_t a, int i) {
  return a->item[i];
}

int darray_index_of(darray_ptr a, void *p);
void darray_remove(darray_ptr a, void *p);
void darray_remove_last(darray_ptr a);
void darray_remove_with_test(darray_ptr a, int (*test)(void *));

/*@manual darray
Removes the pointer at index 'i' in the dynamic array 'a'.
*/
void darray_remove_index(darray_ptr a, int n);
void darray_copy(darray_ptr dst, darray_ptr src);
void darray_remove_all(darray_ptr d);
void darray_forall(darray_t a, void (*func)(void *));
void darray_forall2(darray_t a,
                    void (*func)(void *darray_item, void *scope_ptr),
                    void *scope_ptr);
void darray_forall3(darray_t a,
                    void (*func)(void *darray_item,
                                 void *scope_ptr1,
                                 void *scope_ptr2),
                    void *scope_ptr1,
                    void *scope_ptr2);
void darray_forall4(darray_t a,
                    void (*func)(void *darray_item,
                                 void *scope_ptr1,
                                 void *scope_ptr2,
                                 void *scope_ptr3),
                    void *scope_ptr1,
                    void *scope_ptr2,
                    void *scope_ptr3);

void *darray_at_test(darray_ptr a, int (*test)(void *,void *), void *scope_ptr);

/*@manual darray
Returns the number of pointers held in 'a'.
*/
static inline int darray_count(darray_ptr a) {
  return a->count;
}

static inline int darray_is_empty(darray_ptr a) {
  return !a->count;
}

static inline void *darray_last(darray_t a) {
  return a->item[a->count - 1];
}

#pragma GCC visibility pop

#endif //__PBC_DARRAY_H__
