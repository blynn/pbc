CLEANFILES = *~
maintainer-clean-local:
	-rm -rf Makefile.in

AM_CPPFLAGS = -I../include
LDADD = ../libpbc.la

noinst_PROGRAMS = benchmark testhilbert testbbs testsig testibe \
	testpairing genaparam gendparam geneparam genfparam \
	genbgn listmnt testpoly 19 sing testbls testibs testbb

benchmark_SOURCES = benchmark.c
testhilbert_SOURCES = testhilbert.c
testbbs_SOURCES = testbbs.c
testsig_SOURCES = testsig.c
testibe_SOURCES = testibe.c
testpairing_SOURCES = testpairing.c
genaparam_SOURCES = genaparam.c
gendparam_SOURCES = gendparam.c
geneparam_SOURCES = geneparam.c
genfparam_SOURCES = genfparam.c
genbgn_SOURCES = genbgn.c
listmnt_SOURCES = listmnt.c
testpoly_SOURCES = testpoly.c
19_SOURCES = 19.c
sing_SOURCES = sing.c

testbls_SOURCES = testbls.c
testbls_LDADD = ../libpbc.la ../sig/libpbc_sig.la

testibs_SOURCES = testibs.c
testibs_LDADD = ../libpbc.la ../sig/libpbc_sig.la

testbb_SOURCES = testbb.c
testbb_LDADD = ../libpbc.la ../sig/libpbc_sig.la
