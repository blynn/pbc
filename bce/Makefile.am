CLEANFILES = *~
maintainer-clean-local:
	-rm -rf Makefile.in

AM_CPPFLAGS = -I../include
LDADD = ../libpbc.la

noinst_PROGRAMS = testbce
testbce_SOURCES = testbce.c bce.c
