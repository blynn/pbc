#include <stdlib.h>
#include <stdint.h> // for intptr_t
#include <stdio.h>
#include <string.h>
#include "pbc_utils.h"
#include "pbc_memory.h"

#ifdef SAFE_CLEAN
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

/* pbc_mem is a continuous memory keeping track of its size */
static inline size_t pbc_mem_get_size(size_t *p) {
  return *p;
}

static inline void pbc_mem_set_size(size_t *p, size_t size) {
  *p = size;
}

static inline void *pbc_mem_to_ptr(size_t *p) {
  return p + 1;
}

static inline void *pbc_ptr_to_mem(size_t *p) {
  return p - 1;
}

static void *pbc_mem_malloc(size_t size) {
  void *ptr = malloc(size + sizeof(size_t));
  if(ptr)
    pbc_mem_set_size(ptr, size);
  return ptr;
}

static void pbc_mem_free(void *ptr) {
  memset(ptr, 0, pbc_mem_get_size(ptr) + sizeof(size_t));
  free(ptr);
}

static void *default_pbc_malloc(size_t size) {
  void *ptr = pbc_mem_malloc(size);
  if(!ptr) pbc_die("malloc() error");
  return pbc_mem_to_ptr(ptr);
}

static void *default_pbc_realloc(void *old, size_t new_size) {
  void *new = pbc_mem_malloc(new_size);
  if(!new) pbc_die("realloc() error");
  if(old) {
    old = pbc_ptr_to_mem(old);
    memcpy(pbc_mem_to_ptr(new), pbc_mem_to_ptr(old), pbc_mem_get_size(old));
    pbc_mem_free(old);
  }
  return pbc_mem_to_ptr(new);
}

static void default_pbc_free(void *ptr) {
  if(ptr)
    pbc_mem_free(pbc_ptr_to_mem(ptr));
}
#else
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
#endif

/* release memory got from pbc_malloc only by pbc_free(), do not use free() */
void *(*pbc_malloc)(size_t) = default_pbc_malloc;
/* pbc_realloc guarantees zeroing out the memory before moving old memory */
void *(*pbc_realloc)(void *, size_t) = default_pbc_realloc;
/* pbc_free guarantees zeroing out the memory */
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
