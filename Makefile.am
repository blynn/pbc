SUBDIRS = . example gen

SO_VERSION = 1:0:0

ChangeLog :
	git log > ChangeLog

CLEANFILES = *~
maintainer-clean-local:
	-rm -rf Makefile.in depcomp ltmain.sh config.h.in configure compile config.guess config.sub \
	missing install-sh autom4te.cache aclocal.m4 tmp

lib_LTLIBRARIES = libpbc.la 

# this should really be a versioned dir, i.e., $(includedir)/pbc-0.5.0
library_includedir = $(includedir)/pbc
library_include_HEADERS = include/pbc_a1_param.h \
include/pbc_a_param.h \
include/pbc_curve.h \
include/pbc_d_param.h \
include/pbc_e_param.h \
include/pbc_field.h \
include/pbc_multiz.h \
include/pbc_z.h \
include/pbc_fieldquadratic.h \
include/pbc_f_param.h \
include/pbc_g_param.h \
include/pbc_fp.h \
include/pbc.h \
include/pbc_hilbert.h \
include/pbc_memory.h \
include/pbc_mnt.h \
include/pbc_pairing.h \
include/pbc_param.h \
include/pbc_poly.h \
include/pbc_random.h \
include/pbc_singular.h \
include/pbc_test.h \
include/pbc_utils.h

libpbc_la_CPPFLAGS = -Iinclude
libpbc_la_SOURCES = arith/field.c arith/z.c \
	arith/naivefp.c arith/fastfp.c \
	arith/fp.c arith/fasterfp.c arith/montfp.c \
	arith/multiz.c \
	arith/dlog.c \
	arith/fieldquadratic.c arith/poly.c \
	arith/random.c arith/init_random.c \
	misc/darray.c misc/symtab.c misc/get_time.c \
	misc/utils.c misc/memory.c misc/extend_printf.c \
	ecc/mpc.c ecc/mnt.c ecc/hilbert.c ecc/curve.c ecc/pairing.c \
	ecc/singular.c \
	ecc/param.c ecc/a_param.c ecc/d_param.c ecc/e_param.c \
	ecc/f_param.c ecc/g_param.c
libpbc_la_LDFLAGS = -lgmp -lm -version-info $(SO_VERSION) $(PBC_LDFLAGS) $(LIBPBC_LDFLAGS)

noinst_PROGRAMS = pbc/pbc benchmark/benchmark benchmark/timersa benchmark/ellnet
pbc_pbc_CPPFLAGS = -I include
pbc_pbc_LDADD = libpbc.la
pbc_pbc_SOURCES = pbc/parser.tab.c pbc/lex.yy.c pbc/pbc.c pbc/pbc_getline.c misc/darray.c misc/symtab.c
benchmark_benchmark_CPPFLAGS = -I include
benchmark_benchmark_LDADD = libpbc.la
benchmark_benchmark_SOURCES = benchmark/benchmark.c
benchmark_timersa_CPPFLAGS = -I include
benchmark_timersa_LDADD = libpbc.la
benchmark_timersa_SOURCES = benchmark/timersa.c
benchmark_ellnet_CPPFLAGS = -I include
benchmark_ellnet_LDADD = libpbc.la
benchmark_ellnet_SOURCES = benchmark/ellnet.c
