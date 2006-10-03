#include <assert.h>
#include <stdio.h>
#include <gmp.h>
#include "symtab.h"
#include "fops.h"
#include "parse.h"
#include "tracker.h"
#include "utils.h"
#include <stdlib.h>
#include <string.h>

void param_out_type(FILE *stream, char *s)
{
    fprintf(stream, "type %s\n", s);
}

void param_out_mpz(FILE *stream, char *s, mpz_t z)
{
    fprintf(stream, "%s ", s);
    mpz_out_str(stream, 0, z);
    fprintf(stream, "\n");
}

void param_out_int(FILE *stream, char *s, int i)
{
    mpz_t z;
    mpz_init(z);

    mpz_set_si(z, i);
    param_out_mpz(stream, s, z);
    mpz_clear(z);
}

void param_read_generic (symtab_t tab, fetch_ops_t fops, void *ctx)
{
    assert (fops);
    assert (ctx);
    token_t tok;
    char *s, *s1;

    token_init(tok);
    for (;;) {
	token_get_generic (tok, fops, ctx);
	if (tok->type != token_word) {
	    break;
	}
	s = strdup(tok->s);
	token_get_generic (tok, fops, ctx);
	if (tok->type != token_word) {
	    break;
	}
	s1 = strdup(tok->s);
	symtab_put(tab, s1, s);
	free(s);
    }
    token_clear(tok);
}

void param_read_buf (symtab_t tab, const char *buf, size_t len)
{
    assert (buf);
    tracker_t t;
    tracker_init (&t, buf, len);
    param_read_generic (tab, &fops_buf, &t);
}

void param_read_str (symtab_t tab, FILE *stream) 
{
    assert (stream);
    param_read_generic (tab, &fops_str, stream);
}


void param_clear_tab(symtab_t tab)
{
    UNUSED_VAR (tab);
    //TODO: free table data
}

void lookup_mpz(mpz_t z, symtab_t tab, char *key)
{
    if (!symtab_has(tab, key)) {
	fprintf(stderr, "missing param: `%s'\n", key);
	return;
    }

    char *data = symtab_at(tab, key);

    mpz_set_str(z, data, 0);
}

int lookup_int(symtab_t tab, char *key)
{
    int res;
    mpz_t z;
    if (!symtab_has(tab, key)) {
	fprintf(stderr, "missing param: `%s'\n", key);
	return 0;
    }

    char *data = symtab_at(tab, key);

    mpz_init(z);

    mpz_set_str(z, data, 0);
    res = mpz_get_si(z);

    mpz_clear(z);
    return res;
}
