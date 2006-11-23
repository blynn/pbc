#ifndef __PBC_FOPS_H__
#define __PBC_FOPS_H__

struct fetch_ops_s {
   int (*fops_getc) (void *ctx);
   int (*fops_ungetc) (int c, void *ctx);
};
typedef struct fetch_ops_s fetch_ops_t[1];
typedef struct fetch_ops_s* fetch_ops_ptr;

extern struct fetch_ops_s fops_buf;
extern struct fetch_ops_s fops_str;

#endif
