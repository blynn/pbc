//requires
// * stddef.h
#ifndef __PBC_TRACKER_H__
#define __PBC_TRACKER_H__

typedef struct tracker_t {
    const char *base;
    const char *p;
    const char *end;
} tracker_t;

void tracker_init (tracker_t *t, const char *buf, size_t len);

#endif //__PBC_TRACKER_H__
