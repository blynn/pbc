// requires
// * stdlib.h
#ifndef PBC_MEMORY_H
#define PBC_MEMORY_H
//memory allocation functions
extern void *(*pbc_malloc)(size_t);
extern void *(*pbc_realloc)(void *, size_t);
extern void (*pbc_free)(void *);

void *pbc_calloc(size_t, size_t);

/*@manual alloc
Set custom allocation functions.
The parameters must be function pointers to drop-in replacements for
malloc, realloc and free, except that malloc and realloc should
terminate the program on failure: they must not return.
*/
void pbc_set_memory_functions(void *(*malloc_fn)(size_t),
	void *(*realloc_fn)(void *, size_t), void (*free_fn)(void *));
#endif //PBC_MEMORY_H
