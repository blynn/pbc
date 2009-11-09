#!/bin/bash

# this is all from my ebuild (no doi)

function die() {
    echo mv_things.sh ERROR
    exit 1
}

if [ -z "$1" ]; then
    echo "no dest dir given?"
    exit 1
fi

D=$1
Ex=${D}/usr/share/doc/libpbc0/examples/
Ox=`pwd`
mkdir -p ${Ex}/src

echo "installing examples to ${Ex}"

install -o 0 -g 0 -m 755 gen/genalldparams ${Ex} || die
install -o 0 -g 0 -m 755 benchmark/report_times  ${Ex}/run_tests || die

cp -r param/ ${Ex}/ || die
cp {pbc,benchmark,gen,example}/*.c ${Ex}/src || die
rm ${Ex}/src/*.readline.c || die

find ${Ex} -type d -exec chmod 755 {} \; || die
find ${Ex} -type f -exec chmod 644 {} \; || die

#install -o 0 -g 0 -m 644 exmakefile ${Ex}/src/Makefile
install -o 0 -g 0 -m 755 -d ${D}/usr/bin/ || die

echo "building a real pbc"
(cd pbc; gcc -c pbc_getline.readline.c)
gcc -o realpbc -I. -Iinclude pbc/pbc.c -L .libs -lpbc pbc/pbc_getline.readline.o  -lreadline  pbc_pbc-symtab.o pbc_pbc-parser.tab.o pbc_pbc-darray.o pbc_pbc-lex.yy.o

echo "installing the pbc binary"
install -o 0 -g 0 -m 755 realpbc ${D}/usr/bin/pbc || die
rm realpbc

DEV=`echo ${D}-dev | sed s/libpbc0-dev/libpbc-dev/`
mkdir -p ${DEV}/usr/share/doc/libpbc0
mkdir -p ${DEV}/usr/include
mkdir -p ${DEV}/usr/bin
mkdir -p ${DEV}/usr/lib

mv ${D}/usr/lib/*                      ${DEV}/usr/lib
mv ${DEV}/usr/lib/libpbc*.so.*         ${D}/usr/lib
mv ${D}/usr/include/pbc                ${DEV}/usr/include/
mv ${D}/usr/bin/pbc                    ${DEV}/usr/bin/
mv ${D}/usr/share/doc/libpbc0/examples ${DEV}/usr/share/doc
