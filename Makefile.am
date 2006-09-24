SUBDIRS = . sig bce test 
EXTRA_DIST = Makefile.mingw Makefile.orig param/*.param report_times list_hilbert genalldparams setup \
	bce/bce.h \
	misc/get_time.win32.c \
	doc/pretty/pretty.css doc/makeover doc/extract doc/*.x[ms]l doc/default.css doc/Makefile

ChangeLog :
	cg-log > ChangeLog

CLEANFILES = *~
maintainer-clean-local:
	-rm -rf Makefile.in depcomp ltmain.sh config.h.in configure compile config.guess config.sub \
		missing install-sh autom4te.cache aclocal.m4 tmp key2.stor system.stor

lib_LTLIBRARIES = libpbc.la 

# this should really be a versioned dir, i.e., $(includedir)/pbc-0.2.16/pbc
library_includedir = $(includedir)/pbc
library_include_HEADERS = include/a_param.h include/d_param.h include/darray.h include/f_param.h \
	include/fieldquadratic.h include/hash.h include/leak.h include/mpc.h include/pairing.h \
	include/parse.h include/poly.h include/a1_param.h include/curve.h include/e_param.h \
	include/field.h include/get_time.h include/hilbert.h include/mnt.h include/noleak.h \
	include/param.h include/pbc.h include/random.h include/symtab.h include/utils.h \
	include/fops.h include/tracker.h include/strclone.h

libpbc_la_CPPFLAGS = -Iinclude
libpbc_la_SOURCES = arith/field.c arith/naivefp.c arith/fieldquadratic.c arith/poly.c \
	arith/random.c misc/parse.c misc/darray.c misc/symtab.c misc/get_time.c \
	misc/fops.c ecc/mpc.c ecc/mnt.c ecc/hilbert.c ecc/curve.c ecc/pairing.c \
	ecc/param.c ecc/a1_param.c ecc/a_param.c ecc/d_param.c ecc/e_param.c \
	ecc/f_param.c misc/tracker.c misc/extend_printf.c
libpbc_la_LDFLAGS =  -lgmp -lm

