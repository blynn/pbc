#include <stdlib.h>
#include <stdint.h> // for intptr_t
#include <stdio.h>
#include <string.h>
#include "pbc_utils.h"
#include "pbc_memory.h"

/* guarantee zeroing the memory */
static void gmp_free(void *ptr, size_t size) {
  if(ptr)
    memset(ptr, 0, size);
  free(ptr);
}

static void* gmp_malloc(size_t size) {
  return malloc(size);
}

/* guarantee zeroing the memory
 * realloc() is not suitable for use with secure memory
 * because memory contents are not zeroed out. */
static void* gmp_realloc(void *old_ptr, size_t old_size, size_t new_size) {
  void *new_ptr = malloc(new_size);
  if(new_ptr && old_ptr)
    memcpy(new_ptr, old_ptr, old_size);
  gmp_free(old_ptr, old_size);
  return new_ptr;
}

static void gmp_guarantee_zero_memory(void) {
  __gmp_set_memory_functions(gmp_malloc, gmp_realloc, gmp_free);
}

__attribute__((constructor)) void init(void) {
  gmp_guarantee_zero_memory();
}

static void *default_pbc_malloc(size_t size) {
  void *res = malloc(size);
  if (!res) pbc_die("malloc() error");
  return res;
}

static void *default_pbc_realloc(void *ptr, size_t size) {
  void *res = realloc(ptr, size);
  if (!res) pbc_die("realloc() error");
  return res;
}

static void default_pbc_free(void *ptr) { free(ptr); }

void *(*pbc_malloc)(size_t) = default_pbc_malloc;
void *(*pbc_realloc)(void *, size_t) = default_pbc_realloc;
void (*pbc_free)(void *) = default_pbc_free;

void pbc_set_memory_functions(void *(*malloc_fn)(size_t),
    void *(*realloc_fn)(void *, size_t), void (*free_fn)(void *)) {
  pbc_malloc = malloc_fn;
  pbc_realloc = realloc_fn;
  pbc_free = free_fn;
}

void *pbc_calloc(size_t nmemb, size_t size) {
  void *res = pbc_malloc(nmemb * size);
  if (!res) pbc_die("calloc() error");
  memset(res, 0, nmemb * size);
  return res;
}

char *pbc_strdup(const char *s) {
  int len = strlen(s);
  char *res = pbc_malloc(len + 1);
  strcpy(res, s);
  return res;
}
