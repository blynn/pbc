CLEANFILES = *~
maintainer-clean-local:
	-rm -rf Makefile.in

AM_CPPFLAGS = -I../include
LDADD = ../libpbc.la

noinst_PROGRAMS = benchmark testhilbert testbbs testsig testibe \
	testpairing genaparam gendparam geneparam genfparam \
	gena1param listmnt testpoly 19 sing testbls testibs testbb \
	checkfp timefp testexp pbc timersa yuanli

yuanli_SOURCES = yuanli.c
pbc_SOURCES = pbc.c
checkfp_SOURCES = checkfp.c
timefp_SOURCES = timefp.c
timersa_SOURCES = timersa.c
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
gena1param_SOURCES = gena1param.c
listmnt_SOURCES = listmnt.c
testpoly_SOURCES = testpoly.c
19_SOURCES = 19.c
sing_SOURCES = sing.c
testexp_SOURCES = testexp.c

testbls_SOURCES = testbls.c
testbls_LDADD = ../libpbc.la ../sig/libpbc_sig.la

testibs_SOURCES = testibs.c
testibs_LDADD = ../libpbc.la ../sig/libpbc_sig.la

testbb_SOURCES = testbb.c
testbb_LDADD = ../libpbc.la ../sig/libpbc_sig.la
