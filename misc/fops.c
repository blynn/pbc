#include <assert.h>
#include <stdio.h>
#include "fops.h"
#include "tracker.h"

static int getc_from_buf (void *ctx);
static int ungetc_into_buf (int c, void *ctx);
static int getc_from_str (void *ctx);
static int ungetc_into_str (int c, void *ctx);

struct fetch_ops_s fops_buf = { getc_from_buf, ungetc_into_buf };
struct fetch_ops_s fops_str = { getc_from_str, ungetc_into_str };

static int
getc_from_buf (void *ctx)
{
   assert (ctx);
   tracker_t *t = (tracker_t *) ctx;
   if (t->p == t->end) { return EOF; }
   return *t->p++;
}

static int
ungetc_into_buf (int c, void *ctx)
{
   assert (ctx);
   tracker_t *t = (tracker_t *) ctx;
   if (c == EOF) { return EOF; }
   if (t->base == t->p) { return EOF; }
   t->p--;
   return c;
}

static int
getc_from_str (void *ctx)
{
   assert (ctx);
   return fgetc ((FILE *) ctx);
}

static int
ungetc_into_str (int c, void *ctx)
{
   assert (ctx);
   return ungetc (c, (FILE *) ctx);
}
