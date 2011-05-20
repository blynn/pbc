// Exercises a bug reported by Michael Adjedj.
//
// In ecc/param.c, token_get() would increment a pointer past a terminating
// NUL, so the parser would keep attempting to read key/value pairs for a
// symbol table. If the memory after the string contains a duplicate key,
// then we have a memory leak because we strdup the value and misc/symtab.c
// overwrites existing elements during insert.
//
// Run with valgrind to spot the bug.
#include "pbc.h"

int main(void) {
  pairing_t p;
  pairing_init_set_str(p,
"type a\n"
"q 8780710799663312522437781984754049815806883199414208211028653399266475630880222957078625179422662221423155858769582317459277713367317481324925129998224791\n"
"h 12016012264891146079388821366740534204802954401251311822919615131047207289359704531102844802183906537786776\n"
"r 730750818665451621361119245571504901405976559617\n"
"exp2 159\n"
"exp1 107\n"
"sign1 1\n"
"sign0 1\0a b a b\n"
  );
  pairing_clear(p);
  return 0;
}
