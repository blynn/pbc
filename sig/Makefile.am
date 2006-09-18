CLEANFILES = *~
maintainer-clean-local:
	-rm -rf Makefile.in

AM_CPPFLAGS = -I../include

lib_LTLIBRARIES = libpbc_sig.la
include_HEADERS = ../include/pbc_sig.h

# this should really be a versioned dir, i.e., $(includedir)/pbc-0.2.16/pbc
library_includedir = $(includedir)/pbc
library_include_HEADERS = ../include/pbc_sig.h

libpbc_sig_la_SOURCES = hash.c sig.c
