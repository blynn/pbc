# I use this Makefile rather than the autotools for simplicity and speed.
# Of course, it is less portable. Much of it is due to Hovav Shacham.

.PHONY: target binaries test clean

target: out libpbc.a binaries

CC := gcc
RANLIB := ranlib
warnflags := -Wall -W -Wfloat-equal -Wendif-labels -Wshadow \
	     -Wpointer-arith -Wcast-align -Wstrict-prototypes \
             -Wredundant-decls #-std=c99 -pedantic
CPPFLAGS := -Iinclude
optflags := -O3 -g -pipe -ffast-math -fomit-frame-pointer
LDLIBS := -lgmp -lm
CFLAGS := $(optflags) $(warnflags)

ifeq ($(PLATFORM),win32)
  nonlinux := .win32
  exe_suffix := .exe
  CC := i586-mingw32msvc-gcc
  AR := i586-mingw32msvc-ar
  RANLIB := i586-mingw32msvc-ranlib
  CPPFLAGS := $(CPPFLAGS) -I/home/blynn/cross/gmp/include
  LDFLAGS := -L/home/blynn/cross/gmp/lib
else
  # tcmalloc is faster than normal malloc.
  LDLIBS := $(LDLIBS) -ltcmalloc
endif

libpbc_srcs := \
  $(addsuffix .c,$(addprefix arith/, \
    field fp montfp naivefp fastfp fasterfp fieldmpz fieldquadratic poly \
    random dlog)) \
  $(addsuffix .c,$(addprefix ecc/, \
    curve singular pairing param_utils \
    a_param d_param e_param f_param g_param \
    hilbert mnt mpc)) \
  $(addsuffix .c,$(addprefix misc/, \
    pbc_assert \
    darray symtab \
    parse \
    fops tracker extend_printf memory)) \
  $(addsuffix $(nonlinux).c,misc/get_time arith/init_random)

libpbc_objs := $(libpbc_srcs:.c=.o)

bin_srcs := \
  $(addsuffix .c,$(addprefix example/, \
    bls hess joux paterson yuanli zhangkim zss)) \
  $(addsuffix .c,$(addprefix gen/, \
    gena1param genaparam gendparam geneparam genfparam \
    hilbertpoly listmnt listfreeman)) \
  benchmark/benchmark.c benchmark/timersa.c benchmark/ellnet.c

define demo_tmpl
  examples += out/$(basename $(notdir $(1)))$(exe_suffix)
  out/$(basename $(notdir $(1)))$(exe_suffix) : $(1) libpbc.a ; \
    $(CC) -o $$@ $(LDFLAGS) $$^ $(LOADLIBES) $(LDLIBS)
endef

$(foreach x,$(bin_srcs:.c=.o),$(eval $(call demo_tmpl,$(x))))

pbc_objs := pbc/pbc.o pbc/pbc_getline.readline.o

pbc/pbc : $(pbc_objs) libpbc.a
	$(CC) -o $@ $(LDFLAGS) $^ $(LOADLIBES) $(LDLIBS) -lreadline

binaries : $(examples) pbc/pbc

test_srcs := \
  $(addsuffix .c,$(addprefix guru/, \
    fp_test quadratic_test poly_test testexp testpairing))

tests := $(test_srcs:.c=)

# Object files needed to test Fp.
fp_objs := $(addsuffix .o, \
  arith/field arith/fp arith/naivefp arith/fastfp arith/fasterfp arith/montfp arith/random arith/init_random misc/extend_printf misc/memory misc/pbc_assert)

guru/testexp: guru/testexp.o libpbc.a
guru/testpairing: guru/testpairing.o libpbc.a
guru/fp_test: guru/fp_test.o $(fp_objs)
guru/poly_test: guru/poly_test.o $(fp_objs) arith/poly.o misc/darray.o
guru/quadratic_test: guru/quadratic_test.o $(fp_objs) arith/fieldquadratic.o

test : $(tests)

out: ; mkdir out

srcs := $(libpbc_srcs) $(bin_srcs) $(test_srcs)
objs := $(srcs:.c=.o) $(pbc_objs)

clean: ; -rm -r out $(objs) libpbc.a

# File dependencies for library-making.
# See GNU Make manual, sect. 11.2.
libpbc.a: libpbc.a($(libpbc_objs))
	$(RANLIB) $@
