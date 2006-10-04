// Pairing-Based Calculator
// mainly for demonstration purposes
//
// It's times like these I wish C had garbage collection

#include <string.h>
#include <ctype.h>
#include "pbc.h"
#include "fieldmpz.h"
#include "utils.h"

/* It's much nicer with readline
#include <readline/readline.h>
#include <readline/history.h>
*/

enum {
    t_none = 0,
    t_id,
    t_int,
    t_comma,
    t_lparen,
    t_rparen,
    t_add,
    t_sub,
    t_mul,
    t_div,
    t_set,
    t_pow,
    t_unk,
    t_function,
    t_pairing,
    t_element,
    t_field,
};

enum {
    pe_expect_factor = 100,
    pe_expect_rparen,
    pe_arglist,
    re_varnotfound = 200,
    re_badlvalue,
    re_funnotfound,
};

static field_t Z;

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

void tree_delete(tree_ptr t)
{
    void delete_child(void *p) { tree_delete(p); }
    darray_forall(t->child, delete_child);
    darray_clear(t->child);
    switch(t->type) {
	case t_id:
	case t_function:
	case t_int:
	    id_delete(t->data);
	    break;
    }
    free(t);
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

    if (isdigit(c)) {
	tok_type = t_int;
	word[0] = c;

	int i = 1;
	for (;;) {
	    c = *lexcp++;
	    if (isdigit(c)) {
		word[i++] = c;
	    } else {
		word[i] = '\0';
		lexcp--;
		break;
	    }
	}
	return;
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
	case ',':
	    tok_type = t_comma;
	    break;
	case '=':
	    tok_type = t_set;
	    break;
	case '^':
	    tok_type = t_pow;
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

static int lastparseerror;
static void setparseerror(int i)
{
    lastparseerror = i;
}

static int lastruntimeerror;
static void setruntimeerror(int i)
{
    lastruntimeerror = i;
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
    if (!c) return NULL;
    darray_append(t->child, c);
    for (;;) {
	if (tok_type == t_rparen) {
	    lex();
	    return t;
	}
	if (tok_type != t_comma) {
	    setparseerror(pe_arglist);
	    return NULL;
	}
	lex(); //expect comma
	c = parsesetexpr();
	if (!c) return NULL;
	darray_append(t->child, c);
    }
}

static tree_ptr parsesubfactor(void)
{
    tree_ptr t;
    switch(tok_type) {
	id_ptr id;
	case t_id:
	    id = id_new(word);
	    lex();
	    if (tok_type == t_lparen) {
		if (parseexprlist(t = tree_new(t_function, id))) {
		    return t;
		}
		tree_delete(t);
		return NULL;
	    } else {
		return tree_new(t_id, id);
	    }
	case t_sub:
	    lex();
	    t = parsesubfactor();
	    if (!t) return NULL;
	    tree_ptr t1 = tree_new(t_function, id_new("neg"));
	    darray_append(t1->child, t);
	    return t1;
	case t_lparen:
	    lex();
	    t = parsesetexpr();
	    if (!t) return NULL;
	    if (tok_type != t_rparen) {
		tree_delete(t);
		setparseerror(pe_expect_rparen);
		return NULL;
	    }
	    lex();
	    return t;
	case t_int:
	    id = id_new(word);
	    lex();
	    return tree_new(t_int, id);
	default:
	    setparseerror(pe_expect_factor);
	    return NULL;
    }
}

static tree_ptr parsefactor(void)
{
    tree_ptr t1;
    t1 = parsesubfactor();
    if (!t1) return NULL;
    
    if (tok_type == t_pow) {
	tree_ptr t2, res;
	lex();
	t2 = parsefactor();
	if (!t2) {
	    tree_delete(t1);
	    return NULL;
	}
	res = tree_new(t_function, id_new("pow"));
	darray_append(res->child, t1);
	darray_append(res->child, t2);
	return res;
    }
    return t1;
}

static tree_ptr parseterm(void)
{
    tree_ptr t1, t2, res;
    res = parsefactor();
    if (!res) return NULL;
    for (;;) {
	switch(tok_type) {
	    case t_mul:
		lex();
		t2 = parsefactor();
		if (!t2) {
		    tree_delete(res);
		    return NULL;
		}
		//t1 = tree_new(t_mul, NULL);
		t1 = tree_new(t_function, id_new("mul"));
		darray_append(t1->child, res);
		darray_append(t1->child, t2);
		res = t1;
		break;
	    case t_div:
		lex();
		t2 = parsefactor();
		if (!t2) {
		    tree_delete(res);
		    return NULL;
		}
		//t1 = tree_new(t_div, NULL);
		t1 = tree_new(t_function, id_new("div"));
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
    if (!res) {
	return NULL;
    }
    for (;;) {
	switch(tok_type) {
	    case t_add:
		lex();
		t2 = parseterm();
		if (!t2) {
		    tree_delete(res);
		    return NULL;
		}
		//t1 = tree_new(t_add, NULL);
		t1 = tree_new(t_function, id_new("add"));
		darray_append(t1->child, res);
		darray_append(t1->child, t2);
		res = t1;
		break;
	    case t_sub:
		lex();
		t2 = parseterm();
		if (!t2) {
		    tree_delete(res);
		    return NULL;
		}
		//t1 = tree_new(t_sub, NULL);
		t1 = tree_new(t_function, id_new("sub"));
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
    if (!t1) return NULL;
    if (tok_type == t_set) {
	lex();
	t2 = parsesetexpr();
	if (!t2) {
	    tree_delete(t1);
	    return NULL;
	}
	res = tree_new(t_set, NULL);
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
	printf("NULL");
	return;
    }
    switch (t->type) {
	case t_set:
	    print_tree(t->child->item[0]);
	    printf(" = ");
	    print_tree(t->child->item[1]);
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
		if (i < t->child->count - 1) printf(", ");
	    }
	    printf(")");
	    break;
	default:
	    printf("?!?");
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
	//current policy: always clear elements, always copy elements
	res->data = malloc(sizeof(element_t));
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

static char *aparam =
"type a\n\
q 8780710799663312522437781984754049815806883199414208211028653399266475630880222957078625179422662221423155858769582317459277713367317481324925129998224791\n\
h 12016012264891146079388821366740534204802954401251311822919615131047207289359704531102844802183906537786776\n\
r 730750818665451621361119245571504901405976559617\n\
exp2 159\n\
exp1 107\n\
sign1 1\n\
sign0 1\n";
static val_ptr f_pairing_new_a_default(darray_ptr arg)
{
    UNUSED_VAR(arg);
    pairing_ptr p = malloc(sizeof(pairing_t));
    pairing_init_inp_buf(p, aparam, strlen(aparam));
    return val_new(t_pairing, p);
}

static val_ptr f_pairing_get_group(
	field_ptr (*get_group)(pairing_ptr p), darray_ptr arg)
{
    val_ptr res;
    if (arg->count != 1) {
	printf("expect one argument\n");
	return NULL;
    }
    val_ptr a0 = arg->item[0];
    if (a0->type != t_pairing) {
	return NULL;
    }
    pairing_ptr pairing = a0->data;
    res = val_new(t_field, get_group(pairing));
    return res;
}

static val_ptr f_pairing_G1(darray_ptr arg)
{
    field_ptr getG1(pairing_ptr p) { return p->G1; }
    return f_pairing_get_group(getG1, arg);
}

static val_ptr f_pairing_G2(darray_ptr arg)
{
    field_ptr getG2(pairing_ptr p) { return p->G2; }
    return f_pairing_get_group(getG2, arg);
}

static val_ptr f_pairing_GT(darray_ptr arg)
{
    field_ptr getGT(pairing_ptr p) { return p->GT; }
    return f_pairing_get_group(getGT, arg);
}

static val_ptr f_pairing_Zr(darray_ptr arg)
{
    field_ptr getZr(pairing_ptr p) { return p->Zr; }
    return f_pairing_get_group(getZr, arg);
}

static val_ptr f_random(darray_ptr arg)
{
    val_ptr res;
    if (arg->count != 1) {
	printf("expect one argument\n");
	return NULL;
    }
    val_ptr a0 = arg->item[0];
    if (a0->type != t_field) {
	printf("arg not field!\n");
	return NULL;
    }
    field_ptr f = a0->data;
    element_ptr e = malloc(sizeof(element_t));
    element_init(e, f);
    element_random(e);
    res = val_new(t_element, e);
    return res;
}

static val_ptr f_unary(
	void (*unary)(element_ptr, element_ptr), darray_ptr arg)
{
    val_ptr res;
    if (arg->count != 1) {
	printf("expect one argument\n");
	return NULL;
    }
    val_ptr a0 = arg->item[0];
    if (a0->type != t_element) {
	printf("arg not element!\n");
	return NULL;
    }
    element_ptr e0 = a0->data;
    element_ptr e = malloc(sizeof(element_t));
    element_init(e, e0->field);
    unary(e, e0);
    res = val_new(t_element, e);
    return res;
}

static val_ptr f_bin_op(
	void (*binop)(element_ptr, element_ptr, element_ptr),
	darray_ptr arg)
{
    val_ptr res;
    if (arg->count != 2) {
	printf("expect two arguments\n");
	return NULL;
    }
    val_ptr a0 = arg->item[0];
    val_ptr a1 = arg->item[1];
    if (a0->type != t_element) {
	printf("left arg not element!\n");
	return NULL;
    }
    if (a1->type != t_element) {
	printf("right arg not element!\n");
	return NULL;
    }
    element_ptr e0 = a0->data;
    element_ptr e1 = a1->data;
    if (e0->field != e1->field) {
	printf("field mismatch!\n");
	return NULL;
    }
    element_ptr e = malloc(sizeof(element_t));
    element_init(e, e0->field);
    binop(e, e0, e1);
    res = val_new(t_element, e);
    return res;
}


static val_ptr f_add(darray_ptr arg)
{
    return f_bin_op(element_add, arg);
}

static val_ptr f_mul(darray_ptr arg)
{
    return f_bin_op(element_mul, arg);
}

static val_ptr f_sub(darray_ptr arg)
{
    return f_bin_op(element_sub, arg);
}

static val_ptr f_div(darray_ptr arg)
{
    void invertandmul(element_ptr c, element_ptr a, element_ptr b)
    {
	element_t tmp;
	element_init(tmp, b->field);
	element_invert(tmp, b);
	element_mul(c, a, tmp);
	element_clear(tmp);
    }

    return f_bin_op(invertandmul, arg);
}

static val_ptr f_inv(darray_ptr arg)
{
    return f_unary(element_invert, arg);
}

static val_ptr f_neg(darray_ptr arg)
{
    return f_unary(element_neg, arg);
}

static val_ptr f_pow(darray_ptr arg)
{
    val_ptr res;
    if (arg->count != 2) {
	printf("expect two arguments\n");
	return NULL;
    }
    val_ptr a0 = arg->item[0];
    val_ptr a1 = arg->item[1];
    if (a0->type != t_element) {
	printf("left arg not element!\n");
	return NULL;
    }
    if (a1->type != t_element) {
	printf("right arg not element!\n");
	return NULL;
    }
    element_ptr e0 = a0->data;
    element_ptr e1 = a1->data;
    element_ptr e = malloc(sizeof(element_t));
    mpz_t z;
    mpz_init(z);
    element_to_mpz(z, e1);
    element_init(e, e0->field);
    element_pow_mpz(e, e0, z);
    res = val_new(t_element, e);
    mpz_clear(z);
    return res;
}
static val_ptr f_pairing(darray_ptr arg)
{
    val_ptr res;
    if (arg->count != 3) {
	printf("expect three arguments\n");
	return NULL;
    }
    val_ptr a0 = arg->item[0];
    val_ptr a1 = arg->item[1];
    val_ptr a2 = arg->item[2];
    if (a0->type != t_element) {
	printf("arg 1 not element!\n");
	return NULL;
    }
    if (a1->type != t_element) {
	printf("arg 2 not element!\n");
	return NULL;
    }
    if (a2->type != t_pairing) {
	printf("arg 3 not pairing!\n");
	return NULL;
    }
    element_ptr e0 = a0->data;
    element_ptr e1 = a1->data;
    pairing_ptr p = a2->data;
    if (e0->field != p->G1) {
	printf("arg 1 not from G1!\n");
	return NULL;
    }
    if (e1->field != p->G2) {
	printf("arg 2 not from G2!\n");
	return NULL;
    }
    element_ptr e = malloc(sizeof(element_t));
    element_init(e, p->GT);
    pairing_apply(e, e0, e1, p);
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
	    v = symtab_at(var, id->data);
	    if (!v) {
		setruntimeerror(re_varnotfound);
		return NULL;
	    }
	    return val_copy(v);
	case t_set:
	    t1 = t->child->item[0];
	    if (t1->type != t_id) {
		setruntimeerror(re_badlvalue);
		return NULL;
	    }
	    t2 = t->child->item[1];
	    v = execute_tree(t2);
	    if (!v) return NULL;
	    id = t1->data;
	    // clear what's there first
	    if ((res = symtab_at(var, id->data))) {
		val_delete(res);
	    }
	    symtab_put(var, v, id->data);
	    v = symtab_at(var, id->data);
	    return val_copy(v);
	case t_function:
	    id = t->data;
	    fn = symtab_at(builtin, id->data);
	    if (!fn) {
		setruntimeerror(re_funnotfound);
		return NULL;
	    }
	    darray_init(arg);
	    for (i=0; i<t->child->count; i++) {
		v = execute_tree(t->child->item[i]);
		if (!v) {
		    darray_forall(arg, (void (*)(void *)) val_delete);
		    return NULL;
		}
		darray_append(arg, v);
	    }
	    res = fn(arg);
	    for (i=0; i<arg->count; i++) {
		val_delete(arg->item[i]);
	    }
	    darray_clear(arg);
	    return res;
	case t_int:
	    id = t->data;
	    char *cp;
	    mpz_t z;
	    mpz_init(z);
	    for (cp = id->data; *cp; cp++) {
		mpz_mul_ui(z, z, 10);
		mpz_add_ui(z, z, *cp - '0');
	    }
	    element_ptr e = malloc(sizeof(element_t));
	    element_init(e, Z);
	    element_set_mpz(e, z);
	    mpz_clear(z);
	    return val_new(t_element, e);
	default:
	    return NULL;
    }
}

static void parseline(char *line)
{
    val_ptr v;

    tree_ptr t;
    lexcp = line;
    lex();
    if (tok_type == t_none) return;
    t = parsesetexpr();
    if (0) {
	print_tree(t);
	printf("\n");
    }
    if (t) {
	v = execute_tree(t);
	if (!v) {
	    printf("runtime error (error code = %d)\n", lastruntimeerror);
	} else {
	    if (t->type != t_set) val_print(v);
	    val_delete(v);
	}
	tree_delete(t);
    } else {
	printf("parse error (error code = %d)\n", lastparseerror);
    }
}

int main(void)
{
    symtab_init(var);
    symtab_init(builtin);

    pairing_ptr p = malloc(sizeof(pairing_t));
    pairing_init_inp_buf(p, aparam, strlen(aparam));
    symtab_put(var, val_new(t_pairing, p), "A");
    symtab_put(var, val_new(t_field, p->G1), "G1");
    symtab_put(var, val_new(t_field, p->G2), "G2");
    symtab_put(var, val_new(t_field, p->GT), "GT");
    symtab_put(var, val_new(t_field, p->Zr), "Zr");

    symtab_put(builtin, f_pairing_new_a_default, "pairing_new_a_default");
    symtab_put(builtin, f_pairing_G1, "get_G1");
    symtab_put(builtin, f_pairing_G2, "get_G2");
    symtab_put(builtin, f_pairing_GT, "get_GT");
    symtab_put(builtin, f_pairing_Zr, "get_Zr");
    symtab_put(builtin, f_random, "random");
    symtab_put(builtin, f_random, "rand");
    symtab_put(builtin, f_random, "rnd");
    symtab_put(builtin, f_neg, "neg");
    symtab_put(builtin, f_sub, "sub");
    symtab_put(builtin, f_add, "add");
    symtab_put(builtin, f_pow, "pow");
    symtab_put(builtin, f_mul, "mul");
    symtab_put(builtin, f_inv, "inv");
    symtab_put(builtin, f_inv, "invert");
    symtab_put(builtin, f_div, "div");
    symtab_put(builtin, f_pairing, "pairing");

    field_init_z(Z);

    fprintf(stderr, "Pairing-Based Calculator\n");

    for (;;) {
	char s[1024];
	fgets(s, 1024, stdin);
	if (feof(stdin)) break;
	parseline(s);
	/* readline version:
	char *line = readline(NULL);
	if (!line) break;
	parseline(line);
	if (*line) add_history(line);
	free(line);
	*/
    }
    return 0;
}
