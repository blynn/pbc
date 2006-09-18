#undef malloc
#undef realloc
#undef free
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "darray.h"

struct meminfo_s {
    void *ptr;
    char *file;
    int line;
    size_t size;
};
typedef struct meminfo_s meminfo_t[1];
typedef struct meminfo_s *meminfo_ptr;

static int first = 1;
static darray_t blocklistfn_darray;
static darray_ptr blocklist() {
    if (first) {
	first = 0;
	darray_init(blocklistfn_darray);
    }
    return blocklistfn_darray;
}

void mem_report()
{
    size_t total = 0;
    void print(void *data) {
	meminfo_ptr p = data;
	printf("%s: %d: %d\n", p->file, p->line, p->size);
	total += p->size;
    }

    darray_forall(blocklist(), print);
    printf("total %d\n", total);
}

void *mymalloc(size_t size, char *file, int line)
{
    meminfo_ptr p = malloc(sizeof(meminfo_t));
    p->ptr = malloc(size);
    p->size = size;
    p->file = file;
    p->line = line;
    darray_append(blocklist(), p);
    return p->ptr;
}

void *myrealloc(void *ptr, size_t size, char *file, int line)
{
    int i, n;
    darray_ptr a = blocklist();
    if (!ptr) {
	return mymalloc(size, file, line);
    }
    n = a->count;
    for (i=0; i<n; i++) {
	meminfo_ptr p = a->item[i];
	if (p->ptr == ptr) {
	    p->ptr = realloc(ptr, size);
	    p->size = size;
	    p->file = file;
	    p->line = line;
	    return p->ptr;
	}
    }
    printf("realloc BUG! %X\n", (unsigned) ptr);
    mem_report();
    exit(1);
}

void myfree(void *ptr)
{
    int i, n;
    darray_ptr a = blocklist();
    if (!ptr) return;
    n = a->count;
    for (i=0; i<n; i++) {
	meminfo_ptr p = a->item[i];
	if (p->ptr == ptr) {
	    free(p);
	    darray_remove_index(a, i);
	    free(ptr);
	    return;
	}
    }
    printf("free BUG! %X\n", (unsigned) ptr);
    mem_report();
    exit(1);
}

void *wrapmalloc(size_t size)
{
    return mymalloc(size, __FILE__, __LINE__);
}

void *wraprealloc(void *ptr, size_t oldsize, size_t newsize)
{
    return myrealloc(ptr, newsize, __FILE__, __LINE__);
}

void wrapfree(void *ptr, size_t size)
{
    myfree(ptr);
}

void gmp_leak_check(void)
{
    mp_set_memory_functions(wrapmalloc, wraprealloc, wrapfree);
}

void tag_element_init(void *ptr, char *file, int line)
{
    meminfo_ptr p = malloc(sizeof(meminfo_t));
    p->ptr = ptr;
    p->size = 0;
    p->file = file;
    p->line = line;
    darray_append(blocklist(), p);
}

void tag_element_clear(void *ptr)
{
    int i, n;
    darray_ptr a = blocklist();
    if (!ptr) return;
    n = a->count;
    for (i=n-1; i>=0; i--) {
	meminfo_ptr p = a->item[i];
	if (p->ptr == ptr) {
	    darray_remove_index(a, i);
	    return;
	}
    }
    printf("element_clear BUG! %X\n", (unsigned) ptr);
    mem_report();
    exit(1);
}
