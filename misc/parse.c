#include <assert.h>
#include <stdio.h>  // for EOF
#include <string.h> // for strchr
#include <stdlib.h> // for pbc_malloc, pbc_free
#include "pbc_parse.h"
#include "pbc_memory.h"

// TODO: Replace with Flex.
const char *token_get_generic(token_t tok, const char *input) {
  char *buf;
  int n = 32;
  int i;
  int c;

  // Skip whitespace and comments.
  for(;;) {
    for (;;) {
      c = *input++;
      if (c == EOF) {
        tok->type = token_eof;
        return input;
      }
      if (!strchr(" \t\r\n", c)) break;
    }

    if (c == '#') {
      for(;;) {
        c = *input++;
        if (c == EOF) {
          tok->type = token_eof;
          return input;
        }
        if (c == '\n') break;
      }
    } else break;
  }

  if (c == '<') {
    c = *input++;
    if (c == '/') {
      tok->type = token_langleslash;
      return input;
    }
    input--;
    tok->type = token_langle;
    return input;
  } else if (c == '>') {
    tok->type = token_rangle;
    return input;
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
       c = *input++;
       if (c == EOF || strchr(" \t\r\n</>", c)) break;
    }
    buf[i] = 0;
    input--;
    tok->s = buf;
  }
  return input;
}

void token_init(token_t tok) {
   tok->type = token_none;
   tok->s = NULL;
}

void token_clear(token_t tok) {
   pbc_free(tok->s);
}
