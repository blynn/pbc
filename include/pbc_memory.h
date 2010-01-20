// Requires:
// * stdlib.h
#ifndef __PBC_MEMORY_H__
#define __PBC_MEMORY_H__

// Memory allocation functions used by PBC.
extern void *(*pbc_malloc)(size_t);
extern void *(*pbc_realloc)(void *, size_t);
extern void (*pbc_free)(void *);

void *pbc_calloc(size_t, size_t);

/*@manual alloc
Set custom allocation functions.  The parameters must be function pointers to
drop-in replacements for malloc, realloc and free, except that malloc and
realloc should terminate the program on failure: they must not return in this
case.
*/
void pbc_set_memory_functions(void *(*malloc_fn)(size_t),
        void *(*realloc_fn)(void *, size_t), void (*free_fn)(void *));

char *pbc_strdup(const char *s);

#endif //__PBC_MEMORY_H__
