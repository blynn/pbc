%{
#include <stdarg.h>
#include <stdio.h>
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_field.h"

#include "pbc_tree.h"
#define YYSTYPE tree_ptr
#include "parser.tab.h"
%}

%option noyywrap

%%
[0-9]+                  yylval = tree_new_z(yytext); return NUM;
[a-zA-Z_][a-zA-Z0-9_]*  return ID;
:=                      return ASSIGN;
\+                      return PLUS;
-                       return MINUS;
\/                      return DIVIDE;
\*                      return TIMES;
\^                      return POW;
;                       return SEMI;
,                       return SEMI;
\(                      return LPAR;
\)                      return RPAR;
\[                      return LSQU;
\]                      return RSQU;
[ \t\r\n]*  // Skip whitespace.
.                       return UNKNOWN;
%%
