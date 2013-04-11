CLEANFILES = *~
maintainer-clean-local:
	-rm -rf Makefile.in

AM_CPPFLAGS = -I../include
LDADD = ../libpbc.la -lgmp

noinst_PROGRAMS = bls hess joux paterson yuanli zhangkim zss

bls_SOURCES = bls.c
hess_SOURCES = hess.c
joux_SOURCES = joux.c
paterson_SOURCES = paterson.c
yuanli_SOURCES = yuanli.c
zhangkim_SOURCES = zhangkim.c
zss_SOURCES = zss.c
