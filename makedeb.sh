#!/bin/bash

set -e

# pbc/lex.yy.c
# pbc/lex.yy.h
# pbc/parser.tab.c
# pbc/parser.tab.h

# bison -d -b pbc/parser $^
# pbc/lex.yy.c : pbc/parser.lex
#     flex -o $@ --header-file=pbc/lex.yy.h $^

if [ ! -f pbc/parser.tab.c -o pbc/parser.y -nt pbc/parser.tab.c ]; then
    bison -d -b pbc/parser pbc/parser.y
fi

if [ ! -f pbc/lex.yy.c -o pbc/parser.lex -nt pbc/lex.yy.c ]; then
    flex -o pbc/lex.yy.c --header-file=pbc/lex.yy.h pbc/parser.lex
fi

exit 0

dpkg-buildpackage -rfakeroot
