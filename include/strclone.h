//requires
// * stdlib.h
// * string.h

static inline char *strclone(const char *src)
{
    char *dst = malloc(strlen(src) + 1);
    if (dst) strcpy(dst, src);
    return dst;
}
