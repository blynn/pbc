#include <stdio.h>
#include <string.h>

#include "pbc_memory.h"

char *pbc_getline(const char *prompt) {
  char s[1024];
  if (prompt) fputs(prompt, stdout);
  if (!fgets(s, 1024, stdin)) return NULL;
  if (feof(stdin)) return NULL;
  /* use strdup rather than pbc_strdup. because
   * 1. readline version of this function uses malloc.
   * 2. pbc_malloc called by pbc_strdup may differ from malloc.
   * here we keep consistency.
   */
  return strdup(s);
}
