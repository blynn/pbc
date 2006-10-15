#include <stdio.h>
#include <string.h>

char *getline(void)
{
    char s[1024];
    fgets(s, 1024, stdin);
    if (feof(stdin)) return NULL;
    return strdup(s);
}
