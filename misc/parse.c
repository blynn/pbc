#include <assert.h>
#include <stdio.h> //for EOF
#include <string.h> //strchr
#include <stdlib.h> //for pbc_malloc, pbc_free
#include "pbc_fops.h"
#include "pbc_parse.h"
#include "pbc_tracker.h"
#include "pbc_memory.h"

void token_get_generic(token_t tok, fetch_ops_t fops, void *ctx) {
  assert(fops);
  assert(ctx);
  char *buf;
  int n = 32;
  int i;
  int c;

skipwhitespace:
  for (;;) {
  c = fops->fops_getc (ctx);
  if (c == EOF) {
     tok->type = token_eof;
     return;
  }
  if (!strchr(" \t\r\n", c)) break;
  }

  if (c == '#') {
  for(;;) {
     c = fops->fops_getc (ctx);
     if (c == EOF) {
     tok->type = token_eof;
     return;
     }
     if (c == '\n') break;
  }
  goto skipwhitespace;
  }

  if (c == '<') {
  c = fops->fops_getc (ctx);
  if (c == '/') {
     tok->type = token_langleslash;
     return;
  }
  fops->fops_ungetc (c, ctx);
  tok->type = token_langle;
  return;
  } else if (c == '>') {
  tok->type = token_rangle;
  return;
  } else {
  tok->type = token_word;
  pbc_free(tok->s);
  buf = (char *) pbc_malloc(n);
  i = 0;
  for (;;) {
     buf[i] = c;
     i++;
     if (i == n) {
     n += 32;
     buf = (char *) pbc_realloc(buf, n);
     }
     c = fops->fops_getc(ctx);
     if (c == EOF || strchr(" \t\r\n</>", c)) break;
  }
  buf[i] = 0;
  fops->fops_ungetc(c, ctx);
  tok->s = buf;
  }
}

void token_init(token_t tok) {
   tok->type = token_none;
   tok->s = NULL;
}

void token_clear(token_t tok) {
   pbc_free(tok->s);
}
