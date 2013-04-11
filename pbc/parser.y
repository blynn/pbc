%{
#include <stdarg.h>
#include <stdio.h>
#include <stdint.h> // for intptr_t
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_field.h"

#include "pbc_tree.h"
#define YYSTYPE tree_ptr
void yyerror(const char *s);
int yylex(void);

#define YY_NO_INPUT
#define YY_NO_UNPUT

extern int option_easy;
%}

%error-verbose
%token DEFINE
%token TERMINATOR
%token NUM ID
%token LPAR RPAR LSQU RSQU LBRACE RBRACE COMMA
%right QUESTION COLON
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
  | input stmt { tree_eval_stmt($2); }
  ;

stmt
  : expr TERMINATOR
  | DEFINE ID LPAR parms RPAR LBRACE stmtlist RBRACE {
      $$ = tree_new_define($2, $4, $7);
    }
  ;

stmtlist
  : { $$ = tree_new_empty_stmt_list(); }  // Empty.
  | stmtlist stmt { tree_append($1, $2); }
  ;

parms
  : { $$ = tree_new_empty_parms(); }  // Empty.
  | parms1
  ;

parms1
  : ID { $$ = tree_new_empty_parms(); tree_append($$, $1); }
  | parms1 COMMA ID { tree_append($1, $3); }
  ;

expr
  : multinomial
  | ID ASSIGN expr   { $$ = tree_new_assign($1, $3); }
  | expr QUESTION expr COLON expr  { $$ = tree_new_ternary($1, $3, $5); }
  | molecule
  | molecule LSQU expr RSQU  { $$ = tree_new_item($1, $3); }
  | expr EQ expr     { $$ = tree_new_eq($1, $3); }
  | expr NE expr     { $$ = tree_new_ne($1, $3); }
  | expr LE expr     { $$ = tree_new_le($1, $3); }
  | expr GE expr     { $$ = tree_new_ge($1, $3); }
  | expr LT expr     { $$ = tree_new_lt($1, $3); }
  | expr T_GT expr   { $$ = tree_new_gt($1, $3); }
  | expr PLUS expr   { $$ = tree_new_add($1, $3); }
  | expr MINUS expr  { $$ = tree_new_sub($1, $3); }
  | expr TIMES expr  { $$ = tree_new_mul($1, $3); }
  | expr DIVIDE expr { $$ = tree_new_div($1, $3); }
  | expr POW expr    { $$ = tree_new_pow($1, $3); }
  | MINUS expr %prec UMINUS  { $$ = tree_new_neg($2); }
  ;

// Not quite atoms.
molecule
  : molecule LPAR exprlist RPAR  { $$ = $3; tree_set_fun($$, $1); }
  | LPAR expr RPAR               { $$ = $2; }
  | ID
  ;

exprlist
  : { $$ = tree_new_funcall(); }  // Empty.
  | nonemptyexprlist
  ;

nonemptyexprlist
  : expr   { tree_append($$ = tree_new_funcall(), $1); }
  | nonemptyexprlist COMMA expr { tree_append($1, $3); }
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
  | sequence COMMA expr { tree_append($1, $3); }
  ;
%%
