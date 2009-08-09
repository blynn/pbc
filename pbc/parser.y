%{
#include <stdarg.h>
#include <stdio.h>
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_field.h"

#include "pbc_tree.h"
#define YYSTYPE tree_ptr
extern int option_easy;
%}

%error-verbose
%token NUM ID
%token LPAR RPAR LSQU RSQU COMMA
%token TERMINATOR
%left EQ NE LT T_GT LE GE
%right ASSIGN
%left PLUS MINUS
%left DIVIDE TIMES
%right UMINUS
%right POW
%token UNKNOWN
%token END 0 "end of file"
%%
input
  : // Empty.
  | input stmt TERMINATOR { tree_eval_stmt($2); }
  ;

stmt
  : { $$ = NULL; }  // Empty.
  | expr
  | assign_expr
  ;

assign_expr
  : ID ASSIGN expr         { $$ = tree_new_assign($1, $3); }
  | ID ASSIGN assign_expr  { $$ = tree_new_assign($1, $3); }
  ;

expr
  : multinomial
  | subexpr LPAR exprlist RPAR  { $$ = $3; tree_set_fun($$, $1); }
  | expr EQ expr     { $$ = tree_new_bin(fun_eq, $1, $3); }
  | expr NE expr     { $$ = tree_new_bin(fun_ne, $1, $3); }
  | expr LE expr     { $$ = tree_new_bin(fun_le, $1, $3); }
  | expr GE expr     { $$ = tree_new_bin(fun_ge, $1, $3); }
  | expr LT expr     { $$ = tree_new_bin(fun_lt, $1, $3); }
  | expr T_GT expr   { $$ = tree_new_bin(fun_gt, $1, $3); }
  | expr PLUS expr   { $$ = tree_new_bin(fun_add, $1, $3); }
  | expr MINUS expr  { $$ = tree_new_bin(fun_sub, $1, $3); }
  | expr TIMES expr  { $$ = tree_new_bin(fun_mul, $1, $3); }
  | expr DIVIDE expr { $$ = tree_new_bin(fun_div, $1, $3); }
  | expr POW expr    { $$ = tree_new_bin(fun_pow, $1, $3); }
  | MINUS expr %prec UMINUS  { $$ = tree_new_uminus($2); }
  | subexpr
  ;

// Expressions that are also valid for function calls.
subexpr
  : LPAR expr RPAR   { $$ = $2 }
  | ID
  ;

exprlist
  : { $$ = tree_new_funcall(); }  // Empty.
  | nonemptyexprlist
  ;

nonemptyexprlist
  : expr  { tree_fun_append($$ = tree_new_funcall(), $1); }
  | nonemptyexprlist COMMA expr  { tree_fun_append($1, $3); }
  ;

multinomial
  : NUM
  | numlist
  ;

numlist
  : LSQU sequence RSQU { $$ = $2; }
  ;

sequence
  : expr { $$ = tree_new_list($1); }
  | sequence COMMA expr { tree_append_multiz($1, $3); }
  ;
%%
