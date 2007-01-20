// Pairing-Based Calculator
// mainly for demonstration purposes
//
// It's times like these I wish C had garbage collection

#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <unistd.h> //for getopt
#include "pbc.h"
#include "pbc_fieldmpz.h"
#include "pbc_fp.h"
#include "pbc_utils.h"

#include "config.h"

char *pbc_getline(void);

enum {
    t_none = 0,
    t_id,
    t_int,
    t_string,
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
    t_err,
};

enum {
    pe_expect_factor = 100,
    pe_expect_rparen,
    pe_arglist,
    re_varnotfound = 200,
    re_badlvalue,
    re_funnotfound,
    re_unimplemented,
    re_badargcount,
    re_badarg,
    re_fieldmismatch,
};

static int option_echo = 0;

static field_t Z;

static int tok_type;
//TODO: dynamic allocation:
static char word[1024];

struct id_s {
    char *data;
    int alloc;
};
typedef struct id_s *id_ptr;

id_ptr id_new(char *id)
{
    id_ptr res = pbc_malloc(sizeof(struct id_s));
    res->alloc = strlen(id) + 1;
    res->data = pbc_malloc(res->alloc);
    strcpy(res->data, id);
    return res;
}

void id_delete(id_ptr id)
{
    pbc_free(id->data);
    pbc_free(id);
}

struct tree_s {
    int type;
    void *data;
    darray_t child;
};
typedef struct tree_s *tree_ptr;

tree_ptr tree_new(int type, void *data)
{
    tree_ptr res = pbc_malloc(sizeof(struct tree_s));
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
	case t_string:
	case t_function:
	case t_int:
	    id_delete(t->data);
	    break;
    }
    pbc_free(t);
}

static char *currentline;
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

    //comments start with '#' and end at a newline
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

    //strings
    if (c == '"') {
	tok_type = t_string;
	int i = 0;
	for (;;) {
	    c = *lexcp++;
	    if (!c) {
		//string continues on next line
		word[i++] = '\n';
		pbc_free(currentline);
		currentline = pbc_getline();
		if (!currentline) break;
		if (option_echo) puts(currentline);
		lexcp = currentline;
		c = *lexcp++;
	    }
	    if (c == '"') {
		break;
	    }
	    word[i++] = c;
	}
	word[i] = '\0';
	return;
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

static tree_ptr parseprimitive(void)
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
	case t_string:
	    lex();
	    return tree_new(t_string, id_new(word));
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
    
static tree_ptr parsepow(void)
{
    tree_ptr t1;
    t1 = parseprimitive();
    if (tok_type == t_pow) {
	tree_ptr t2, res;
	lex();
	t2 = parseprimitive();
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

static tree_ptr parsefactor(void)
{
    tree_ptr t;
    if (tok_type == t_sub) {
	lex();
	t = parsefactor();
	if (!t) return NULL;
	tree_ptr t1 = tree_new(t_function, id_new("neg"));
	darray_append(t1->child, t);
	return t1;
    }

    t = parsepow();
    return t;
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

static int lastruntimeerror;
static val_ptr newruntimeerror(int i)
{
    val_ptr res = pbc_malloc(sizeof(struct val_s));
    lastruntimeerror = i;
    res->type = t_err;
    res->data = int_to_voidp(i);
    return res;
}

val_ptr val_new(int type, void *data)
{
    val_ptr res = pbc_malloc(sizeof(struct val_s));
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
	    field_out_info(stdout, field);
	    break;
	case t_string:
	    printf("%s", (char *) v->data);
	    break;
	default:
	    printf("val type %d unknown\n", v->type);
	    break;
    }
}

val_ptr val_copy(val_ptr v)
{
    val_ptr res = pbc_malloc(sizeof(struct val_s));
    res->type = v->type;
    if (v->type == t_element) {
	//current policy: always clear elements, always copy elements
	res->data = pbc_malloc(sizeof(element_t));
	element_ptr e = v->data;
	element_init(res->data, e->field);
	element_set(res->data, e);
    } else if (v->type == t_string) {
	res->data = pbc_strdup(v->data);
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
	    pbc_free(v->data);
	    break;
	case t_string:
	    pbc_free(v->data);
	    break;
	case t_err:
	    break;
	case t_pairing:
	    break;
	case t_field:
	    break;
	default:
	    printf("val_delete: case %d not handled: memory leak\n", v->type);
	    break;
    }
    pbc_free(v);
}

struct fun_s {
    val_ptr (*f)(darray_ptr);
    int arity;
    int type[32]; //TODO: replace with darray? who needs more than 32 args?
};

typedef val_ptr (*fun)(darray_ptr);

static val_ptr check_arg(darray_ptr arg, int n, ...)
{
    va_list ap;
    int i;
    val_ptr res = NULL;

    va_start(ap, n);
    if (arg->count != n) {
	printf("expect %d argument(s)\n", n);
	res = newruntimeerror(re_badargcount);
    } else for (i=0; i<n; i++) {
	int t = va_arg(ap, int);
	val_ptr vp = arg->item[i];
	if (vp->type != t) {
	    printf("arg not type %d\n", t);
	    return newruntimeerror(re_badarg);
	    break;
	}
    }

    va_end(ap);
    return res;
}

static val_ptr f_pairing_get_group(
	field_ptr (*get_group)(pairing_ptr p), darray_ptr arg)
{
    val_ptr res;
    res = check_arg(arg, 1, t_pairing);
    if (res) return res;
    val_ptr a0 = arg->item[0];
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
    res = check_arg(arg, 1, t_field);
    if (res) return res;
    val_ptr a0 = arg->item[0];
    field_ptr f = a0->data;
    element_ptr e = pbc_malloc(sizeof(element_t));
    element_init(e, f);
    element_random(e);
    res = val_new(t_element, e);
    return res;
}

static val_ptr f_order(darray_ptr arg)
{
    val_ptr res;
    res = check_arg(arg, 1, t_field);
    if (res) return res;
    val_ptr a0 = arg->item[0];
    field_ptr f = a0->data;

    element_ptr e = pbc_malloc(sizeof(element_t));
    element_init(e, Z);
    element_set_mpz(e, f->order);
    res = val_new(t_element, e);
    return res;
}

static val_ptr f_unary(
	void (*unary)(element_ptr, element_ptr), darray_ptr arg)
{
    val_ptr res;
    res = check_arg(arg, 1, t_element);
    if (res) return res;
    val_ptr a0 = arg->item[0];
    element_ptr e0 = a0->data;
    element_ptr e = pbc_malloc(sizeof(element_t));
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
    res = check_arg(arg, 2, t_element, t_element);
    if (res) return res;
    val_ptr a0 = arg->item[0];
    val_ptr a1 = arg->item[1];
    element_ptr e0 = a0->data;
    element_ptr e1 = a1->data;
    if (e0->field != e1->field) {
	printf("field mismatch!\n");
	return newruntimeerror(re_fieldmismatch);
    }
    element_ptr e = pbc_malloc(sizeof(element_t));
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
    return f_bin_op(element_div, arg);
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
    res = check_arg(arg, 2, t_element, t_element);
    if (res) return res;
    val_ptr a0 = arg->item[0];
    val_ptr a1 = arg->item[1];
    element_ptr e0 = a0->data;
    element_ptr e1 = a1->data;
    element_ptr e = pbc_malloc(sizeof(element_t));
    mpz_t z;
    mpz_init(z);
    element_to_mpz(z, e1);
    element_init(e, e0->field);
    element_pow_mpz(e, e0, z);
    res = val_new(t_element, e);
    mpz_clear(z);
    return res;
}

static pairing_ptr current_pairing;
static val_ptr f_pairing(darray_ptr arg)
{
    val_ptr res;
    if (arg->count < 2) {
	printf("expect at least two arguments\n");
	return newruntimeerror(re_badargcount);
    }
    if (arg->count > 3) {
	printf("expect at most three arguments\n");
	return newruntimeerror(re_badargcount);
    }
    val_ptr a0 = arg->item[0];
    val_ptr a1 = arg->item[1];
    val_ptr a2 = arg->item[2];
    if (a0->type != t_element) {
	printf("arg 1 not element!\n");
	return newruntimeerror(re_badarg);
    }
    if (a1->type != t_element) {
	printf("arg 2 not element!\n");
	return newruntimeerror(re_badarg);
    }
    pairing_ptr p;
    if (arg->count == 3) {
	if (a2->type != t_pairing) {
	    printf("arg 3 not pairing!\n");
	    return newruntimeerror(re_badarg);
	} else {
	    p = a2->data;
	}
    } else {
	p = current_pairing;
    }
    element_ptr e0 = a0->data;
    element_ptr e1 = a1->data;
    if (e0->field != p->G1) {
	printf("arg 1 not from G1!\n");
	return newruntimeerror(re_badarg);
    }
    if (e1->field != p->G2) {
	printf("arg 2 not from G2!\n");
	return newruntimeerror(re_badarg);
    }
    element_ptr e = pbc_malloc(sizeof(element_t));
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
		return newruntimeerror(re_varnotfound);
	    }
	    return val_copy(v);
	case t_set:
	    t1 = t->child->item[0];
	    if (t1->type != t_id) {
		return newruntimeerror(re_badlvalue);
	    }
	    t2 = t->child->item[1];
	    v = execute_tree(t2);
	    if (v->type == t_err) return v;
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
		return newruntimeerror(re_funnotfound);
	    }
	    darray_init(arg);
	    for (i=0; i<t->child->count; i++) {
		v = execute_tree(t->child->item[i]);
		if (v->type == t_err) {
		    darray_forall(arg, (void (*)(void *)) val_delete);
		    return v;
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
	    element_ptr e = pbc_malloc(sizeof(element_t));
	    element_init(e, Z);
	    element_set_mpz(e, z);
	    mpz_clear(z);
	    return val_new(t_element, e);
	case t_string:
	    id = t->data;
	    return val_new(t_string, pbc_strdup(id->data));
	default:
	    return newruntimeerror(re_unimplemented);
    }
}

static void parseline(void)
{
    val_ptr v;

    tree_ptr t;
    lex();
    if (tok_type == t_none) return;
    t = parsesetexpr();
    if (0) {
	print_tree(t);
	printf("\n");
    }
    if (t) {
	v = execute_tree(t);
	if (v) {
	    if (v->type == t_err) {
		printf("runtime error (error code = %d)\n", lastruntimeerror);
	    } else {
		if (t->type != t_set) val_print(v);
	    }
	    val_delete(v);
	}
	tree_delete(t);
    } else {
	printf("parse error (error code = %d)\n", lastparseerror);
    }
}

static char *aparam =
"type a\n\
q 8780710799663312522437781984754049815806883199414208211028653399266475630880222957078625179422662221423155858769582317459277713367317481324925129998224791\n\
h 12016012264891146079388821366740534204802954401251311822919615131047207289359704531102844802183906537786776\n\
r 730750818665451621361119245571504901405976559617\n\
exp2 159\n\
exp1 107\n\
sign1 1\n\
sign0 1\n";

static char *dparam =
"type d\n\
q 625852803282871856053922297323874661378036491717\n\
n 625852803282871856053923088432465995634661283063\n\
h 3\n\
r 208617601094290618684641029477488665211553761021\n\
a 581595782028432961150765424293919699975513269268\n\
b 517921465817243828776542439081147840953753552322\n\
k 6\n\
nk 60094290356408407130984161127310078516360031868417968262992864809623507269833854678414046779817844853757026858774966331434198257512457993293271849043664655146443229029069463392046837830267994222789160047337432075266619082657640364986415435746294498140589844832666082434658532589211525696\n\
hk 1380801711862212484403205699005242141541629761433899149236405232528956996854655261075303661691995273080620762287276051361446528504633283152278831183711301329765591450680250000592437612973269056\n\
coeff0 472731500571015189154958232321864199355792223347\n\
coeff1 352243926696145937581894994871017455453604730246\n\
coeff2 289113341693870057212775990719504267185772707305\n\
nqr 431211441436589568382088865288592347194866189652\n";

static char *eparam =
"type e\n\
q 7245986106510086080714203333362098431608853335867425877960916928496629182991629664903654100214900946450053872786629995869445693724001299041657434948257845644905153122838458864000479326695430719258600053239930483226650953770354174712511646273516974069245462534034085895319225452125649979474047163305307830001\n\
r 730750862221594424981965739670091261094297337857\n\
h 13569343110918781839835249021482970252603216587988030044836106948825516930173270978617489032334001006615524543925753725725046733884363846960470444404747241287743773746682188521738728797153760275116924829183670000\n\
a 7130970454025799000067946137594446075551569949583815943390108723282396973737794273397246892274981883807989525599540630855644968426794929215599380425269625872763801485968007136000471718335185787206876242871042697778608875139078711621836858237429403052273312335081163896980825048123655535355411494046493419999\n\
b 7169309004853894693616698536183663527570664411678352588247044791687141043489072737232715961588288238022010974661903752526911876859197052490952065266265699130144252031591491045333807587788600764557450846327338626261289568016170532652061787582791926724597362401398804563093625182790987016728290050466098223333\n\
exp2 159\n\
exp1 135\n\
sign1 1\n\
sign0 1\n";

static char *fparam =
"type f\n\
q 205523667896953300194896352429254920972540065223\n\
r 205523667896953300194895899082072403858390252929\n\
b 40218105156867728698573668525883168222119515413\n\
beta 115334401956802802075595682801335644058796914268\n\
alpha0 191079354656274778837764015557338301375963168470\n\
alpha1 71445317903696340296199556072836940741717506375\n";

static pairing_t pairing_A, pairing_D, pairing_E, pairing_F;

static void set_pairing_groups(pairing_ptr p) {
    symtab_put(var, val_new(t_field, p->G1), "G1");
    symtab_put(var, val_new(t_field, p->G2), "G2");
    symtab_put(var, val_new(t_field, p->GT), "GT");
    symtab_put(var, val_new(t_field, p->Zr), "Zr");
    symtab_put(var, val_new(t_pairing, p), "current_pairing");
    current_pairing = p;
}

static val_ptr f_init_pairing(darray_ptr arg)
{
    val_ptr res;

    res = check_arg(arg, 1, t_pairing);
    if (res) return res;

    val_ptr a0 = arg->item[0];
    pairing_ptr p = a0->data;
    set_pairing_groups(p);
    return NULL;
}

static val_ptr f_nextprime(darray_ptr arg)
{
    mpz_t p;
    val_ptr res;

    res = check_arg(arg, 1, t_element);
    if (res) return res;
    val_ptr a0 = arg->item[0];
    element_ptr e0 = a0->data;
    if (e0->field != Z) {
	printf("arg not integer!\n");
	return newruntimeerror(re_badarg);
    }
    element_ptr e = pbc_malloc(sizeof(element_t));
    element_init(e, Z);
    mpz_init(p);
    element_to_mpz(p, e0);
    mpz_nextprime(p, p);
    element_set_mpz(e, p);
    res = val_new(t_element, e);
    mpz_clear(p);
    return res;
}

static val_ptr f_brute_force_dlog(darray_ptr arg)
{
    val_ptr res;
    res = check_arg(arg, 2, t_element, t_element);
    if (res) return res;
    val_ptr a0 = arg->item[0];
    val_ptr a1 = arg->item[1];
    element_ptr e0 = a0->data;
    element_ptr e1 = a1->data;
    if (e0->field != e1->field) {
	printf("arg field mismatch!\n");
	return newruntimeerror(re_badarg);
    }
    element_ptr e = pbc_malloc(sizeof(element_t));
    element_init(e, Z);
    brute_force_dlog(e, e0, e1);
    res = val_new(t_element, e);
    return res;
}
static val_ptr f_pollard_rho(darray_ptr arg)
{
    val_ptr res;
    res = check_arg(arg, 3, t_element, t_element, t_field);
    if (res) return res;
    val_ptr a0 = arg->item[0];
    val_ptr a1 = arg->item[1];
    val_ptr a2 = arg->item[2];
    element_ptr e0 = a0->data;
    element_ptr e1 = a1->data;
    if (e0->field != e1->field) {
	printf("arg field mismatch!\n");
	return newruntimeerror(re_badarg);
    }
    field_ptr f = a2->data;
    element_ptr e = pbc_malloc(sizeof(element_t));
    element_init(e, f);
    pollard_rho(e, e0, e1);
    res = val_new(t_element, e);
    return res;
}

static val_ptr f_zz(darray_ptr arg)
{
    mpz_t p;
    val_ptr res;
    res = check_arg(arg, 1, t_element);
    if (res) return res;
    val_ptr a0 = arg->item[0];
    element_ptr e0 = a0->data;
    if (e0->field != Z) {
	printf("arg not integer!\n");
	return newruntimeerror(re_badarg);
    }
    field_ptr f = pbc_malloc(sizeof(field_t));
    mpz_init(p);
    element_to_mpz(p, e0);
    field_init_fp(f, p);
    res = val_new(t_field, f);
    mpz_clear(p);
    return res;
}

static val_ptr f_gen_A(darray_ptr arg)
{
    mpz_t rbits, qbits;
    a_param_t param;
    pairing_ptr p;
    val_ptr res;
    res = check_arg(arg, 2, t_element, t_element);
    if (res) return res;
    val_ptr a0 = arg->item[0];
    val_ptr a1 = arg->item[1];
    element_ptr e0 = a0->data;
    if (e0->field != Z) {
	printf("arg not integer!\n");
	return newruntimeerror(re_badarg);
    }
    element_ptr e1 = a1->data;
    if (e1->field != Z) {
	printf("arg not integer!\n");
	return newruntimeerror(re_badarg);
    }
    mpz_init(rbits);
    mpz_init(qbits);
    element_to_mpz(rbits, e0);
    element_to_mpz(qbits, e1);
    //TODO: check rbits and qbits aren't too big
    a_param_init(param);
    a_param_gen(param, mpz_get_ui(rbits), mpz_get_ui(qbits));
    p = pbc_malloc(sizeof(pairing_t));
    pairing_init_a_param(p, param);
    res = val_new(t_pairing, p);
    mpz_clear(rbits);
    mpz_clear(qbits);
    a_param_clear(param);
    return res;
}

static val_ptr f_fromZZ(darray_ptr arg)
{
    val_ptr res;
    res = check_arg(arg, 2, t_element, t_field);
    if (res) return res;
    val_ptr a0 = arg->item[0];
    val_ptr a1 = arg->item[1];
    element_ptr e = a0->data;
    field_ptr f = a1->data;
    if (e->field != Z) {
	printf("arg not integer!\n");
	return newruntimeerror(re_badarg);
    }
    element_ptr e1 = pbc_malloc(sizeof(element_t));
    element_init(e1, f);
    element_set_mpz(e1, e->data);
    res = val_new(t_element, e1);
    return res;
}

static val_ptr f_fromstr(darray_ptr arg)
{
    val_ptr res;
    res = check_arg(arg, 2, t_string, t_field);
    if (res) return res;
    val_ptr a0 = arg->item[0];
    val_ptr a1 = arg->item[1];
    field_ptr f = a1->data;
    element_ptr e1 = pbc_malloc(sizeof(element_t));
    element_init(e1, f);
    element_set_str(e1, a0->data, 0);
    res = val_new(t_element, e1);
    return res;
}

static val_ptr f_index_calculus(darray_ptr arg)
{
    val_ptr res;
    res = check_arg(arg, 2, t_element, t_element);
    if (res) return res;
    val_ptr a0 = arg->item[0];
    val_ptr a1 = arg->item[1];
    element_ptr e0 = a0->data;
    element_ptr e1 = a1->data;
    element_ptr e = pbc_malloc(sizeof(element_t));
    mpz_t x, g, h, q1;

    //TODO: check e0, e1 are from an integer mod ring
    mpz_init(x);
    mpz_init(g);
    mpz_init(h);
    mpz_init(q1);

    mpz_sub_ui(q1, e0->field->order, 1);

    element_init(e, Z);
    element_to_mpz(g, e0);
    element_to_mpz(h, e1);
    index_calculus_dlog(x, g, h, q1);
    element_set_mpz(e, x);
    res = val_new(t_element, e);
    mpz_clear(x);
    mpz_clear(g);
    mpz_clear(h);
    mpz_clear(q1);
    return res;
}

int main(int argc, char **argv)
{
    for (;;)
    {
	int c = getopt(argc, argv, "e");
	if (c == -1) break;
	switch (c) {
	    case 'e':
		option_echo = 1;
		break;
	    default:
		fprintf(stderr, "unrecognized option: %c\n", c);
		break;
	}
    }

    symtab_init(var);
    symtab_init(builtin);

    pairing_init_inp_buf(pairing_A, aparam, strlen(aparam));
    pairing_init_inp_buf(pairing_D, dparam, strlen(dparam));
    pairing_init_inp_buf(pairing_E, eparam, strlen(eparam));
    pairing_init_inp_buf(pairing_F, fparam, strlen(fparam));
    symtab_put(var, val_new(t_pairing, pairing_A), "A");
    symtab_put(var, val_new(t_pairing, pairing_D), "D");
    symtab_put(var, val_new(t_pairing, pairing_E), "E");
    symtab_put(var, val_new(t_pairing, pairing_F), "F");

    set_pairing_groups(pairing_A);

    symtab_put(builtin, f_init_pairing, "init_pairing");
    symtab_put(builtin, f_pairing_G1, "get_G1");
    symtab_put(builtin, f_pairing_G2, "get_G2");
    symtab_put(builtin, f_pairing_GT, "get_GT");
    symtab_put(builtin, f_pairing_Zr, "get_Zr");
    symtab_put(builtin, f_random, "random");
    symtab_put(builtin, f_random, "rand");
    symtab_put(builtin, f_random, "rnd");
    symtab_put(builtin, f_order, "order");
    symtab_put(builtin, f_order, "ord");
    symtab_put(builtin, f_neg, "neg");
    symtab_put(builtin, f_sub, "sub");
    symtab_put(builtin, f_add, "add");
    symtab_put(builtin, f_pow, "pow");
    symtab_put(builtin, f_mul, "mul");
    symtab_put(builtin, f_inv, "inv");
    symtab_put(builtin, f_inv, "invert");
    symtab_put(builtin, f_div, "div");
    symtab_put(builtin, f_pairing, "pairing");
    symtab_put(builtin, f_nextprime, "nextprime");
    symtab_put(builtin, f_brute_force_dlog, "brute_force_dlog");
    symtab_put(builtin, f_pollard_rho, "pollard_rho");
    symtab_put(builtin, f_index_calculus, "index_calculus");
    symtab_put(builtin, f_zz, "ZZ");
    symtab_put(builtin, f_gen_A, "gen_A");
    symtab_put(builtin, f_fromZZ, "fromZZ");
    symtab_put(builtin, f_fromstr, "fromstr");

    field_init_z(Z);

    fprintf(stderr, "%s\n", PACKAGE_STRING);

    for (;;) {
	currentline = pbc_getline();
	if (!currentline) break;
	if (option_echo) puts(currentline);
	lexcp = currentline;
	parseline();
	pbc_free(currentline);
    }
    return 0;
}
