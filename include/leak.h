#include <stdlib.h>

#define malloc(x) mymalloc(x, __FILE__, __LINE__)
#define realloc(x, y) myrealloc(x, y, __FILE__, __LINE__)
#define free myfree

void *mymalloc(size_t size, char *file, int line);
void *myrealloc(void *ptr, size_t size, char *file, int line);
void myfree(void *ptr);

//For GMP. Use someting like
//mp_set_memory_functions(wrapmalloc, wraprealloc, wrapfree);
void *wrapmalloc(size_t size);
void *wraprealloc(void *ptr, size_t oldsize, size_t newsize);
void wrapfree(void *ptr, size_t size);
void gmp_leak_check(void);
void tag_element_init(void *ptr, char *file, int line);
void tag_element_clear(void *ptr);

#define element_init(x,y) realelement_init(x,y); tag_element_init(x, __FILE__, __LINE__)
#define element_clear(x) realelement_clear(x); tag_element_clear(x)

