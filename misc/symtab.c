#include <stdlib.h>
#include <string.h>
#include "pbc_memory.h"

#include "symtab.h"

struct entry_s {
  char *key;
  void *data;
};
typedef struct entry_s *entry_ptr;
typedef struct entry_s entry_t[1];

void symtab_init(symtab_t t) {
  darray_init(t->list);
}

static void clear(void *data) {
  entry_ptr e = data;
  pbc_free(e->key);
  pbc_free(e);
}

void symtab_clear(symtab_t t) {
  darray_forall(t->list, clear);
  darray_clear(t->list);
}

void symtab_put(symtab_t t, void *data, const char *key) {
  int i, n = t->list->count;
  entry_ptr e;
  for (i=0; i<n; i++) {
    e = t->list->item[i];
    if (!strcmp(e->key, key)) goto doit;
  }
  e = pbc_malloc(sizeof(entry_t));
  e->key = pbc_strdup(key);
  darray_append(t->list, e);
doit:
  e->data = data;
}

int symtab_has(symtab_t t, const char *key) {
  int i, n = t->list->count;
  for (i = 0; i < n; i++) {
    entry_ptr e = t->list->item[i];
    if (!strcmp(e->key, key)) return 1;
  }
  return 0;
}

void *symtab_at(symtab_t t, const char *key) {
  int i, n = t->list->count;
  for (i=0; i<n; i++) {
    entry_ptr e = t->list->item[i];
    if (!strcmp(e->key, key)) return e->data;
  }
  return NULL;
}

void symtab_forall_data(symtab_t t, void (*func)(void *)) {
  int i, n = t->list->count;
  for (i=0; i<n; i++) {
    entry_ptr e = t->list->item[i];
    func(e->data);
  }
}
