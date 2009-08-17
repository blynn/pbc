#include <stdio.h>
#include <string.h>

#include "pbc_memory.h"

char *pbc_getline(const char *prompt) {
  char s[1024];
  if (prompt) fputs(prompt, stdout);
  if (!fgets(s, 1024, stdin)) return NULL;
  if (feof(stdin)) return NULL;
  return pbc_strdup(s);
}
