%{
#include <stdarg.h>
#include <stdio.h>
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_field.h"

#include "pbc_tree.h"
#define YYSTYPE tree_ptr
%}

%error-verbose
%token NUM ID
%token LPAR RPAR LSQU RSQU COMMA
%token TERMINATOR
%right ASSIGN
%left PLUS MINUS
%left DIVIDE TIMES
%right UMINUS
%right POW
%token UNKNOWN
%%
input
  : // Empty.
  | input stmt TERMINATOR { tree_eval_stmt($2); }
  ;

stmt
  : expr
  | assign_expr
  ;

assign_expr
  : ID ASSIGN expr         { $$ = tree_new_assign($1, $3); }
  | ID ASSIGN assign_expr  { $$ = tree_new_assign($1, $3); }
  ;

expr
  : NUM
  | ID LPAR exprlist RPAR
  | ID
  | expr PLUS expr   { $$ = tree_new_bin(fun_add, $1, $3); }
  | expr MINUS expr  { $$ = tree_new_bin(fun_sub, $1, $3); }
  | expr TIMES expr  { $$ = tree_new_bin(fun_mul, $1, $3); }
  | expr DIVIDE expr { $$ = tree_new_bin(fun_div, $1, $3); }
  | expr POW expr    { $$ = tree_new_bin(fun_pow, $1, $3); }
  | MINUS expr %prec UMINUS
  | LPAR expr RPAR
  ;

exprlist
  : // Empty.
  | nonemptyexprlist
  ;

nonemptyexprlist
  : expr
  | nonemptyexprlist COMMA expr
  ;
%%
