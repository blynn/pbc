#ifndef __PBC_SYMTAB_H__
#define __PBC_SYMTAB_H__

#include "darray.h"

#pragma GCC visibility push(hidden)

struct symtab_s {
    darray_t list;
};
typedef struct symtab_s symtab_t[1];
typedef struct symtab_s *symtab_ptr;

/*@manual symtab
Initialize symbol table 't'. Must be called before 't' is used.
*/
void symtab_init(symtab_t t);

/*@manual symtab
Clears symbol table 't'. Should be called after 't' is no longer needed.
*/
void symtab_clear(symtab_t t);

/*@manual symtab
Puts 'value' at 'key' in 't'.
*/
void symtab_put(symtab_t t, void *value, const char *key);

/*@manual symtab
Returns true if 't' contains key 'key'.
*/
int symtab_has(symtab_t t, const char *key);

/*@manual symtab
Returns pointer at key 'key' in 't'.
*/
void *symtab_at(symtab_t t, const char *key);

void symtab_forall_data(symtab_t t, void (*func)(void *));

#pragma GCC visibility pop

#endif //__PBC_SYMTAB_H__
