#ifndef __UTIL_H__
#define __UTIL_H__

#include <stdlib.h>
#include <string.h>

#define UNUSED_VAR(a) do { /* nothing */ } while (&a == 0)

static inline char *strclone(const char *src)
{
    char *dst = malloc(strlen(src) + 1);
    if (dst) strcpy(dst, src);
    return dst;
}

#endif
