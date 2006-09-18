#ifndef __FOPS_H__
#define __FOPS_H__

typedef struct fetch_ops_t {
   int (*fops_getc) (void *ctx);
   int (*fops_ungetc) (int c, void *ctx);
} fetch_ops_t;

extern fetch_ops_t fops_buf;
extern fetch_ops_t fops_str;

#endif
