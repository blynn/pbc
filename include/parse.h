#ifndef PARSE_H
#define PARSE_H

#include <stdio.h>
#include "fops.h"

enum {
    token_none = 0,
    token_langle,
    token_langleslash,
    token_rangle,
    token_word,
    token_eof,
};

struct token_s {
    int type;
    char *s;
};
typedef struct token_s token_t[1];
typedef struct token_s *token_ptr;

void token_init(token_t tok);
void token_clear(token_t tok);
void token_get(token_t tok, FILE *stream);
void token_get_from_buf (token_t tok, const char *buf, char **p, size_t len);
void token_get_generic (token_t tok, fetch_ops_t fops, void *ctx);

#endif //PARSE_H
