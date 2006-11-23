//requires
// * fops.h
#ifndef __PBC_PARSE_H__
#define __PBC_PARSE_H__

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
void token_get_generic (token_t tok, fetch_ops_t fops, void *ctx);

#endif //__PBC_PARSE_H__
