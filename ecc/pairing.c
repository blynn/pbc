#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_field.h"
#include "pbc_darray.h"
#include "pbc_poly.h"
#include "pbc_curve.h"
#include "pbc_symtab.h"
#include "pbc_pairing.h"
#include "pbc_memory.h"

#include "pbc_param.h"
#include "pbc_a_param.h"
#include "pbc_mnt.h"
#include "pbc_d_param.h"
#include "pbc_e_param.h"
#include "pbc_f_param.h"
#include "pbc_a1_param.h"
#include "pbc_g_param.h"

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

static const char *token_get(token_t tok, const char *input) {
  char *buf;
  int n = 32;
  int i;
  char c;

  // Skip whitespace and comments.
  for(;;) {
    for (;;) {
      c = *input++;
      if (!c) {
        tok->type = token_eof;
        return input;
      }
      if (!strchr(" \t\r\n", c)) break;
    }

    if (c == '#') {
      for(;;) {
        c = *input++;
        if (!c) {
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
       if (!c || strchr(" \t\r\n</>", c)) break;
    }
    buf[i] = 0;
    input--;
    tok->s = buf;
  }
  return input;
}

static void token_init(token_t tok) {
   tok->type = token_none;
   tok->s = NULL;
}

static void token_clear(token_t tok) {
   pbc_free(tok->s);
}

static void read_symtab(symtab_t tab, const char *input) {
  token_t tok;
  char *s, *s1;

  token_init(tok);
  for (;;) {
    input = token_get(tok, input);
    if (tok->type != token_word) {
      break;
    }
    s = pbc_strdup(tok->s);
    input = token_get(tok, input);
    if (tok->type != token_word) {
      pbc_free(s);
      break;
    }
    s1 = pbc_strdup(tok->s);
    symtab_put(tab, s1, s);
    pbc_free(s);
  }
  token_clear(tok);
}

int generic_is_almost_coddh(element_ptr a, element_ptr b,
    element_ptr c, element_ptr d, pairing_t pairing) {
  int res = 0;
  element_t t0, t1;

  element_init(t0, pairing->GT);
  element_init(t1, pairing->GT);
  pairing_apply(t0, a, d, pairing);
  pairing_apply(t1, b, c, pairing);
  if (!element_cmp(t0, t1)) {
    res = 1;
  } else {
    element_mul(t0, t0, t1);
    if (element_is1(t0)) res = 1;
  }
  element_clear(t0);
  element_clear(t1);
  return res;
}

static void phi_warning(element_ptr out, element_ptr in, pairing_ptr pairing) {
  UNUSED_VAR(out);
  UNUSED_VAR(in);
  UNUSED_VAR(pairing);
  printf("Phi() not implemented for this pairing type yet!\n");
}

static void default_option_set(struct pairing_s *pairing, char *key, char *value) {
  UNUSED_VAR(pairing);
  UNUSED_VAR(key);
  UNUSED_VAR(value);
}

static void default_pp_init(pairing_pp_t p, element_ptr in1, pairing_t pairing) {
  UNUSED_VAR(pairing);
  p->data = (void *) in1;
}

static void default_pp_apply(element_ptr out, element_ptr in2, pairing_pp_t p) {
  p->pairing->map(out, p->data, in2, p->pairing);
}

static void default_pp_clear(pairing_pp_t p) {
  UNUSED_VAR(p);
}

int pairing_init_set_str(pairing_t pairing, const char *input) {
  pairing->option_set = default_option_set;
  pairing->pp_init = default_pp_init;
  pairing->pp_clear = default_pp_clear;
  pairing->pp_apply = default_pp_apply;
  pairing->is_almost_coddh = generic_is_almost_coddh;
  pairing->phi = phi_warning;

  symtab_t tab;
  symtab_init(tab);
  read_symtab(tab, input);
  const char *lookup(const char *key) {
    if (!symtab_has(tab, key)) {
      pbc_error("missing param: `%s'", key);
      return NULL;
    }
    return symtab_at(tab, key);
  }
  const char *s = lookup("type");

  static struct {
    char *s;
    void (*fun)(pbc_param_ptr, const char *(*)(const char *));
  } funtab[] = {
      { "a", pbc_param_init_a },
      { "d", pbc_param_init_d },
      { "e", pbc_param_init_e },
      { "f", pbc_param_init_f },
      { "g", pbc_param_init_g },
      { "a1", pbc_param_init_a1 },
  };

  int res = 1;
  if (s) {
    unsigned int i;
    for(i = 0; i < sizeof(funtab)/sizeof(*funtab); i++) {
      if (!strcmp(s, funtab[i].s)) {
	pbc_param_t par;
	funtab[i].fun(par, lookup);
	pairing_init_pbc_param(pairing, par);
	pbc_param_clear(par);
	res = 0;
	break;
      }
    }
  }

  symtab_forall_data(tab, pbc_free);
  symtab_clear(tab);
  if (res) {
    pbc_error("unknown pairing type");
    return 1;
  }

  pairing->G1->pairing = pairing;
  pairing->G2->pairing = pairing;
  pairing->GT->pairing = pairing;
  return 0;
}

void pairing_clear(pairing_t pairing) {
  pairing->clear_func(pairing);
}

// TODO: it's most likely better to add extra stuff to field_t
// so no new data structures are needed to create mulitplicative subgroups.
// Additionally the same code could be used with curve_t
// Will consider it later, especially if timings turn out bad

static void gt_out_info(FILE *out, field_ptr f) {
  gmp_fprintf(out, "GT: order %Zd roots of unity in:\n", f->order);
  field_out_info(out, f->data);
}

static void gt_from_hash(element_ptr e, void *data, int len) {
  pairing_ptr pairing = e->field->pairing;
  element_from_hash(e->data, data, len);
  pairing->finalpow(e);
}

static void gt_random(element_ptr e) {
  pairing_ptr pairing = e->field->pairing;
  element_random(e->data);
  pairing->finalpow(e);
}

// multiplicative subgroup of a field
static void mulg_field_clear(field_t f) {
  UNUSED_VAR(f);
}

static void mulg_init(element_ptr e) {
  e->data = pbc_malloc(sizeof(element_t));
  field_ptr f = e->field->data;
  element_init(e->data, f);
  element_set1(e->data);
}

static void mulg_clear(element_ptr e) {
  element_clear(e->data);
  pbc_free(e->data);
}

static void mulg_set(element_ptr x, element_t a) {
  element_set(x->data, a->data);
}

static int mulg_cmp(element_ptr x, element_t a) {
  return element_cmp(x->data, a->data);
}

static size_t mulg_out_str(FILE *stream, int base, element_ptr e) {
  return element_out_str(stream, base, e->data);
}

static int mulg_set_str(element_ptr e, char *s, int base) {
  return element_set_str(e->data, s, base);
}

static int mulg_to_bytes(unsigned char *data, element_ptr e) {
  return element_to_bytes(data, e->data);
}

static int mulg_from_bytes(element_ptr e, unsigned char *data) {
  return element_from_bytes(e->data, data);
}

static int mulg_length_in_bytes(element_ptr e) {
  return element_length_in_bytes(e->data);
}

static int mulg_snprint(char *s, size_t n, element_ptr e) {
  return element_snprint(s, n, e->data);
}

static void mulg_to_mpz(mpz_ptr z, element_ptr e) {
  element_to_mpz(z, e->data);
}

static void mulg_set1(element_t e) {
  element_set1(e->data);
}

static void mulg_mul(element_ptr x, element_t a, element_t b) {
  element_mul(x->data, a->data, b->data);
}

static void mulg_div(element_ptr x, element_t a, element_t b) {
  element_div(x->data, a->data, b->data);
}

static void mulg_invert(element_ptr x, element_t a) {
  element_invert(x->data, a->data);
}

static int mulg_is1(element_ptr x) {
  return element_is1(x->data);
}

static void mulg_pow_mpz(element_t x, element_t a, mpz_t n) {
  element_pow_mpz(x->data, a->data, n);
}

static void mulg_pp_init(element_pp_t p, element_t in) {
  p->data = pbc_malloc(sizeof(element_pp_t));
  element_pp_init(p->data, in->data);
}

static void mulg_pp_clear(element_pp_t p) {
  element_pp_clear(p->data);
  pbc_free(p->data);
}

static void mulg_pp_pow(element_t out, mpz_ptr power, element_pp_t p) {
  element_pp_pow(out->data, power, p->data);
}

void GT_init_finite_field(pairing_ptr pairing, field_t f) {
  field_ptr gt = pairing->GT;

  field_init(gt);
  gt->data = f;
  f->pairing = pairing;

  mpz_set(gt->order, pairing->r);
  gt->field_clear = mulg_field_clear;
  gt->out_info = gt_out_info;

  gt->init = mulg_init;
  gt->clear = mulg_clear;
  gt->set = mulg_set;
  gt->cmp = mulg_cmp;

  gt->out_str = mulg_out_str;
  gt->set_str = mulg_set_str;
  gt->to_bytes = mulg_to_bytes;
  gt->from_bytes = mulg_from_bytes;
  gt->length_in_bytes = mulg_length_in_bytes;
  gt->fixed_length_in_bytes = f->fixed_length_in_bytes;
  gt->to_mpz = mulg_to_mpz;
  gt->snprint = mulg_snprint;

  // TODO: set gt->nqr to something?
  // set is_sqr, sqrt to something?

  // additive notation
  gt->set0 = mulg_set1;
  gt->add = mulg_mul;
  gt->sub = mulg_div;
  gt->mul_mpz = mulg_pow_mpz;
  gt->neg = mulg_invert;
  gt->is0 = mulg_is1;

  // multiplicative notation
  gt->set1 = mulg_set1;
  gt->mul = mulg_mul;
  gt->div = mulg_div;
  gt->pow_mpz = mulg_pow_mpz;
  gt->invert = mulg_invert;
  gt->is1 = mulg_is1;
  gt->pp_init = mulg_pp_init;
  gt->pp_clear = mulg_pp_clear;
  gt->pp_pow = mulg_pp_pow;

  gt->random = gt_random;
  gt->from_hash = gt_from_hash;
}
