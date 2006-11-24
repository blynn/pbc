CLEANFILES = *~
maintainer-clean-local:
	-rm -rf Makefile.in

AM_CPPFLAGS = -I../include
LDADD = ../libpbc.la

noinst_PROGRAMS = bls bonehfranklin hess joux paterson yuanli zhang zss

bls_SOURCES = bls.c
bonehfranklin_SOURCES = bonehfranklin.c
hess_SOURCES = hess.c
joux_SOURCES = joux.c
paterson_SOURCES = paterson.c
yuanli_SOURCES = yuanli.c
zhang_SOURCES = zhang.c
zss_SOURCES = zss.c
