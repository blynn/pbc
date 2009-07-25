#include <stdio.h>
#include <string.h>

#include "pbc_memory.h"

char *pbc_getline(void)
{
    char s[1024];
    if (!fgets(s, 1024, stdin)) return NULL;
    if (feof(stdin)) return NULL;
    return pbc_strdup(s);
}
