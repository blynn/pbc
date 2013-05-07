#include <stdio.h>
#include <stdlib.h>
#include <stdint.h> // for intptr_t
#include <string.h>
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_memory.h"
#include "pbc_param.h"
#include "pbc_a_param.h"
#include "pbc_mnt.h"
#include "pbc_d_param.h"
#include "pbc_e_param.h"
#include "pbc_f_param.h"
#include "pbc_a1_param.h"
#include "pbc_g_param.h"
#include "pbc_i_param.h"

#include "misc/symtab.h"
#include "ecc/param.h"

// Parser that reads a bunch of strings and places them in a symbol table.
// TODO: Replace with Flex/Bison?

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

// Reads next token from `input`.
// Returns 1 on reaching `end` (if not NULL) or '\0' is read, 0 otherwise.
static const char *token_get(token_t tok, const char *input, const char *end) {
  char *buf;
  int n = 32;
  int i;
  char c;
  #define get() (((!end || input < end) && *input) ? (c = *input++, 0) : 1)
  // Skip whitespace and comments.
  for(;;) {
    do {
      if (get()) {
        tok->type = token_eof;
        return input;
      }
    } while (strchr(" \t\r\n", c));
    if (c == '#') {
      do {
        if (get()) {
          tok->type = token_eof;
          return input;
        }
      } while (c != '\n');
    } else break;
  }

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
    if (get() || strchr(" \t\r\n</>", c)) break;
  }
  buf[i] = 0;
  tok->s = buf;
  return input;
  #undef get
}

static void token_init(token_t tok) {
   tok->type = token_none;
   tok->s = NULL;
}

static void token_clear(token_t tok) {
   pbc_free(tok->s);
}

static void read_symtab(symtab_t tab, const char *input, size_t limit) {
  token_t tok;
  const char *inputend = limit ? input + limit : NULL;
  token_init(tok);
  for (;;) {
    input = token_get(tok, input, inputend);
    if (tok->type != token_word) break;
    char *key = pbc_strdup(tok->s);
    input = token_get(tok, input, inputend);
    if (tok->type != token_word) {
      pbc_free(key);
      break;
    }
    symtab_put(tab, pbc_strdup(tok->s), key);
    pbc_free(key);
  }
  token_clear(tok);
}

// These functions have hidden visibility (see header).

void param_out_type(FILE *stream, char *s) {
  fprintf(stream, "type %s\n", s);
}

void param_out_mpz(FILE *stream, char *s, mpz_t z) {
  fprintf(stream, "%s ", s);
  mpz_out_str(stream, 0, z);
  fprintf(stream, "\n");
}

void param_out_int(FILE *stream, char *s, int i) {
  mpz_t z;
  mpz_init(z);

  mpz_set_si(z, i);
  param_out_mpz(stream, s, z);
  mpz_clear(z);
}

static const char *lookup(symtab_t tab, const char *key) {
  if (!symtab_has(tab, key)) {
    pbc_error("missing param: `%s'", key);
    return NULL;
  }
  return symtab_at(tab, key);
}

int lookup_mpz(mpz_t z, symtab_t tab, const char *key) {
  const char *data = lookup(tab, key);
  if (!data) {
    pbc_error("missing param: `%s'", key);
    return 1;
  }
  mpz_set_str(z, data, 0);
  return 0;
}

int lookup_int(int *n, symtab_t tab, const char *key) {
  mpz_t z;
  const char *data = lookup(tab, key);
  if (!data) {
    pbc_error("missing param: `%s'", key);
    return 1;
  }
  mpz_init(z);

  mpz_set_str(z, data, 0);
  *n = mpz_get_si(z);
  mpz_clear(z);

  return 0;
}

static int param_set_tab(pbc_param_t par, symtab_t tab) {
  const char *s = lookup(tab, "type");

  static struct {
    char *s;
    int (*fun)(pbc_param_ptr, symtab_t tab);
  } funtab[] = {
      { "a", pbc_param_init_a },
      { "d", pbc_param_init_d },
      { "e", pbc_param_init_e },
      { "f", pbc_param_init_f },
      { "g", pbc_param_init_g },
      { "a1", pbc_param_init_a1 },
      { "i", pbc_param_init_i },
  };

  int res = 1;
  if (s) {
    unsigned int i;
    for(i = 0; i < sizeof(funtab)/sizeof(*funtab); i++) {
      if (!strcmp(s, funtab[i].s)) {
        res = funtab[i].fun(par, tab);
        if (res) pbc_error("bad pairing parameters");
        return res;
      }
    }
  }

  pbc_error("unknown pairing type");
  return res;
}

// Public functions:

int pbc_param_init_set_str(pbc_param_t par, const char *input) {
  symtab_t tab;
  symtab_init(tab);
  read_symtab(tab, input, 0);
  int res = param_set_tab(par, tab);
  symtab_forall_data(tab, pbc_free);
  symtab_clear(tab);
  return res;
}

int pbc_param_init_set_buf(pbc_param_t par, const char *input, size_t len) {
  symtab_t tab;
  symtab_init(tab);
  read_symtab(tab, input, len);
  int res = param_set_tab(par, tab);
  symtab_forall_data(tab, pbc_free);
  symtab_clear(tab);
  return res;
}
