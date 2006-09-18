#ifndef __TRACKER_H__
#define __TRACKER_H__

#include <stdlib.h>

typedef struct tracker_t {
    const char *base;
    const char *p;
    const char *end;
} tracker_t;

void tracker_init (tracker_t *t, const char *buf, size_t len);

#endif
