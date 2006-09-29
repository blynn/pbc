// Pairing-Based Calculator
// meant for demos only
// long-term plan: write wrappers for scripting language like Lua
//
// It's times like these I wish C had garbage collection

#include <string.h>
#include <ctype.h>
#include "pbc.h"
#include "utils.h"

enum {
    t_none = 0,
    t_id,
    t_comma,
    t_lparen,
    t_rparen,
    t_add,
    t_sub,
    t_mul,
    t_div,
    t_set,
    t_unk,
    t_function,
    t_pairing,
    t_element,
    t_field,
};

static int tok_type;
static char word[128];

struct id_s {
    char *data;
    int alloc;
};
typedef struct id_s *id_ptr;

id_ptr id_new(char *id)
{
    id_ptr res = malloc(sizeof(struct id_s));
    res->alloc = strlen(id) + 1;
    res->data = malloc(res->alloc);
    strcpy(res->data, id);
    return res;
}

void id_delete(id_ptr id)
{
    free(id->data);
    free(id);
}

struct tree_s {
    int type;
    void *data;
    darray_t child;
};
typedef struct tree_s *tree_ptr;

tree_ptr tree_new(int type, void *data)
{
    tree_ptr res = malloc(sizeof(struct tree_s));
    res->type = type;
    res->data = data;
    darray_init(res->child);
    return res;
}

static char *lexcp;

static void lex(void)
{
    char c;
    if (!lexcp) {
	tok_type = t_none;
	return;
    }
    c = *lexcp++;
  skipwhitespace:
    for (;;) {
	if (!strchr(" \t\r\n", c)) break;
	if (!c) {
	    tok_type = t_none;
	    return;
	}
	c = *lexcp++;
    }

    if (c == '#') {
	for (;;) {
	    c = *lexcp++;
	    if (!c) {
		tok_type = t_none;
		return;
	    }
	    if (c == '\n') break;
	}
	goto skipwhitespace;
    }

    if (isalpha(c) || c == '_') {
	tok_type = t_id;
	word[0] = c;

	int i = 1;
	for (;;) {
	    c = *lexcp++;
	    if (isalnum(c) || c == '_') {
		word[i++] = c;
	    } else {
		word[i] = '\0';
		lexcp--;
		break;
	    }
	}
	return;
    }

    switch(c) {
	case '.':
	    tok_type = t_comma;
	    break;
	case '=':
	    tok_type = t_set;
	    break;
	case '*':
	    tok_type = t_mul;
	    break;
	case '/':
	    tok_type = t_div;
	    break;
	case '+':
	    tok_type = t_add;
	    break;
	case '-':
	    tok_type = t_sub;
	    break;
	case '(':
	    tok_type = t_lparen;
	    break;
	case ')':
	    tok_type = t_rparen;
	    break;
	default:
	    tok_type = t_unk;
	    break;
    }
}

static tree_ptr parsesetexpr(void);

static tree_ptr parseexprlist(tree_ptr t)
{
    tree_ptr c;
    lex(); // expect lparen
    if (tok_type == t_rparen) {
	lex();
	return t;
    }
    c = parsesetexpr();
    darray_append(t->child, c);
    for (;;) {
	if (tok_type == t_rparen) {
	    lex();
	    return t;
	}
	lex(); //expect comma
	c = parsesetexpr();
	darray_append(t->child, c);
    }
}

static tree_ptr parsefactor(void)
{
    tree_ptr t;
    switch(tok_type) {
	id_ptr id;
	case t_id:
	    lex();
	    id = id_new(word);
	    if (tok_type == t_lparen) {
		parseexprlist(t = tree_new(t_function, id));
		return t;
	    } else {
		return tree_new(t_id, id);
	    }
	case t_lparen:
	    lex();
	    t = parsesetexpr();
	    lex(); //expect rparen
	    return t;
	default:
	    return NULL;
    }
}

static tree_ptr parseterm(void)
{
    tree_ptr t1, t2, res;
    res = parsefactor();
    for (;;) {
	switch(tok_type) {
	    case t_mul:
		lex();
		t1 = tree_new(t_mul, NULL);
		t2 = parsefactor();
		darray_append(t1->child, res);
		darray_append(t1->child, t2);
		res = t1;
		break;
	    case t_div:
		lex();
		t1 = tree_new(t_div, NULL);
		t2 = parsefactor();
		darray_append(t1->child, res);
		darray_append(t1->child, t2);
		res = t1;
		break;
	    default:
		return res;
	}
    }
}

static tree_ptr parseexpr(void)
{
    tree_ptr t1, t2, res;
    res = parseterm();
    for (;;) {
	switch(tok_type) {
	    case t_add:
		lex();
		t1 = tree_new(t_add, NULL);
		t2 = parseterm();
		darray_append(t1->child, res);
		darray_append(t1->child, t2);
		res = t1;
		break;
	    case t_sub:
		lex();
		t1 = tree_new(t_sub, NULL);
		t2 = parseterm();
		darray_append(t1->child, res);
		darray_append(t1->child, t2);
		res = t1;
		break;
	    default:
		return res;
	}
    }
}

static tree_ptr parsesetexpr(void)
{
    tree_ptr t1, t2, res;
    t1 = parseexpr();
    if (tok_type == t_set) {
	lex();
	res = tree_new(t_set, NULL);
	t2 = parsesetexpr();
	darray_append(res->child, t1);
	darray_append(res->child, t2);
	return res;
    }
    return t1;
}

static void print_tree(tree_ptr t)
{
    id_ptr id;
    int i;
    if (!t) {
	printf("NULL\n");
	return;
    }
    switch (t->type) {
	case t_set:
	    print_tree(t->child->item[0]);
	    printf(" = ");
	    print_tree(t->child->item[1]);
	    printf("\n");
	    break;
	case t_id:
	    id = t->data;
	    printf("%s", id->data);
	    break;
	case t_function:
	    id = t->data;
	    printf("%s(", id->data);
	    for (i=0; i<t->child->count; i++) {
		print_tree(t->child->item[i]);
	    }
	    printf(")\n");
	    break;
	default:
	    printf("?!?\n");
	    break;
    }
}

static symtab_t var;
static symtab_t builtin;

struct val_s {
    int type;
    void *data;
};
typedef struct val_s *val_ptr;

val_ptr val_new(int type, void *data)
{
    val_ptr res = malloc(sizeof(struct val_s));
    res->type = type;
    res->data = data;
    return res;
}

static void val_print(val_ptr v)
{
    pairing_ptr pairing;
    field_ptr field;
    element_ptr e;
    switch (v->type) {
	case t_element:
	    e = v->data;
	    element_out_str(stdout, 0, e);
	    printf("\n");
	    break;
	case t_pairing:
	    pairing = v->data;
	    printf("pairing: G1bits=%d G2bits=%d GTbits=%d\n",
		    pairing_length_in_bytes_x_only_G1(pairing) * 8,
		    pairing_length_in_bytes_x_only_G2(pairing) * 8,
		    pairing_length_in_bytes_GT(pairing) * 8);
	    break;
	case t_field:
	    field = v->data;
	    field_print_info(stdout, field);
	    break;
	default:
	    printf("val type %d unknown\n", v->type);
	    break;
    }
}

val_ptr val_copy(val_ptr v)
{
    val_ptr res = malloc(sizeof(struct val_s));
    res->type = v->type;
    if (v->type == t_element) {
	element_ptr e = v->data;
	element_init(res->data, e->field);
	element_set(res->data, e);
    } else {
	res->data = v->data;
    }

    return res;
}

void val_delete(val_ptr v)
{
    switch(v->type) {
	case t_element:
	    //current policy: always clear elements, always copy elements
	    element_clear(v->data);
	    free(v->data);
	    break;
	case t_pairing:
	    //TODO: can't do stuff like this until refcounting is implemented.
	    //pairing_clear(v->data);
	    //free(v->data);
	    break;
	case t_field:
	    break;
	default:
	    printf("val_delete: case %d not handled: memory leak\n", v->type);
	    break;
    }
}

struct fun_s {
    val_ptr (*f)(darray_ptr);
    int arity;
    int type[32]; //TODO: replace with darray? who needs more than 32 args?
};

typedef val_ptr (*fun)(darray_ptr);
static val_ptr f_pairing_new_a_default(darray_ptr arg)
{
    UNUSED_VAR(arg);
    static char *param =
"type a\n\
q 8780710799663312522437781984754049815806883199414208211028653399266475630880222957078625179422662221423155858769582317459277713367317481324925129998224791\n\
h 12016012264891146079388821366740534204802954401251311822919615131047207289359704531102844802183906537786776\n\
r 730750818665451621361119245571504901405976559617\n\
exp2 159\n\
exp1 107\n\
sign1 1\n\
sign0 1\n";
    pairing_ptr p = malloc(sizeof(pairing_t));
    pairing_init_inp_buf(p, param, strlen(param));
    return val_new(t_pairing, p);
}

static val_ptr f_pairing_G1(darray_ptr arg)
{
    val_ptr res;
    val_ptr a0 = arg->item[0];
    pairing_ptr pairing = a0->data;
    res = val_new(t_field, pairing->G1);
    return res;
}

static val_ptr f_pairing_G2(darray_ptr arg)
{
    val_ptr res;
    val_ptr a0 = arg->item[0];
    pairing_ptr pairing = a0->data;
    res = val_new(t_field, pairing->G2);
    return res;
}

static val_ptr f_pairing_GT(darray_ptr arg)
{
    val_ptr res;
    val_ptr a0 = arg->item[0];
    pairing_ptr pairing = a0->data;
    res = val_new(t_field, pairing->GT);
    return res;
}

static val_ptr f_random(darray_ptr arg)
{
    val_ptr res;
    val_ptr a0 = arg->item[0];
    field_ptr f = a0->data;
    element_ptr e = malloc(sizeof(element_t));
    element_init(e, f);
    element_random(e);
    res = val_new(t_element, e);
    return res;
}

static val_ptr execute_tree(tree_ptr t)
{
    darray_t arg;
    id_ptr id;
    fun fn;
    int i;
    val_ptr res, v;
    tree_ptr t1, t2;

    switch (t->type) {
	case t_id:
	    id = t->data;
	    printf("var lookup %s\n", id->data);
	    v = symtab_at(var, id->data);
	    printf("var type %d\n", v->type);
	    if (!v) {
		printf("no such variable\n");
		res = NULL;
		break;
	    }
	    res = val_copy(v);
	    printf("var lookup done\n");
	    break;
	case t_set:
	    printf("set begin\n");
	    t1 = t->child->item[0];
	    if (t1->type != t_id) {
		printf("invalid lvalue\n");
		res = NULL;
		break;
	    }
	    t2 = t->child->item[1];
printf("set exec\n");
	    v = execute_tree(t2);
printf("set symput\n");
	    id = t1->data;
	    symtab_put(var, v, id->data);
printf("set copy\n");
	    res = val_copy(v);
printf("set done\n");
	    break;
	case t_function:
	    id = t->data;
printf("function %s\n", id->data);
	    fn = symtab_at(builtin, id->data);
	    if (!fn) {
		printf("no such function %s\n", id->data);
		return NULL;
	    }
	    darray_init(arg);
	    for (i=0; i<t->child->count; i++) {
		darray_append(arg, execute_tree(t->child->item[i]));
	    }
printf("fun call\n");
	    res = fn(arg);
printf("fun ret\n");
	    for (i=0; i<arg->count; i++) {
		val_delete(arg->item[i]);
	    }
	    darray_clear(arg);
	    break;
	default:
	    res = NULL;
	    break;
    }
    return res;
}

static void parseline(char *line)
{
    val_ptr v;

    tree_ptr t;
    lexcp = line;
    lex();
    t = parsesetexpr();
    if (1) print_tree(t);
    v = execute_tree(t);
    if (!v) {
	printf("error\n");
    } else {
	val_print(v);
	val_delete(v);
    }
}

int main(void)
{
    char s[1024];
    symtab_init(var);
    symtab_init(builtin);

    symtab_put(builtin, f_pairing_new_a_default, "pairing_new_a_default");
    symtab_put(builtin, f_pairing_G1, "pairing_G1");
    symtab_put(builtin, f_pairing_G2, "pairing_G2");
    symtab_put(builtin, f_pairing_GT, "pairing_GT");
    symtab_put(builtin, f_random, "random");
    printf("Pairing-Based Calculator\n");

    for (;;) {
	gets(s); //TODO: the dreaded gets(), useful for rapid development though
	if (feof(stdin)) break;
	parseline(s);
    }
    return 0;
}
