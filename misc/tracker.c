#include <assert.h>
#include <stddef.h>
#include "pbc_tracker.h"

void tracker_init (tracker_t *t, const char *buf, size_t len)
{
    assert (t);
    assert (buf);
    t->base = buf;
    t->p = buf;
    t->end = buf+len;
}
