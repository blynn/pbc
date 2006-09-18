#ifndef SYMTAB_H
#define SYMTAB_H

#include "darray.h"

struct symtab_s {
    darray_t list;
};
typedef struct symtab_s symtab_t[1];
typedef struct symtab_s symtab_ptr;

void symtab_init(symtab_t t);
void symtab_clear(symtab_t t);
void symtab_put(symtab_t t, void *data, char *key);
int symtab_has(symtab_t t, char *key);
void *symtab_at(symtab_t t, char *key);

#endif //SYMTAB_H
