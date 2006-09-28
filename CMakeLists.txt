cmake_minimum_required(VERSION 2.4)
project(pbc)
include_directories(./include)
find_library(GMP_LIBRARY NAMES gmp PATHS /usr/lib /usr/local/lib)
set(leaksrc misc/leak.c)
set(unused
foo
)
add_library(pbc
arith/montfp.c
arith/naivefp.c
arith/tinyfp.c arith/fastfp.c arith/fasterfp.c
arith/field.c
arith/fp.c
arith/fieldquadratic.c
arith/poly.c arith/random.c
ecc/curve.c ecc/pairing.c ecc/param.c
ecc/a1_param.c ecc/a_param.c ecc/d_param.c ecc/e_param.c ecc/f_param.c
ecc/hilbert.c ecc/mnt.c ecc/mpc.c
misc/darray.c misc/symtab.c
misc/parse.c
misc/fops.c misc/tracker.c
misc/extend_printf.c
misc/get_time.c
)
add_library(pbc_sig
sig/sig.c sig/hash.c
)
foreach(bin
test/19          test/geneparam  test/testbbs      test/testibe
test/benchmark   test/genfparam  test/testbls      test/testibs
test/gena1param  test/listmnt    test/testfi       test/testpairing
test/genaparam   test/sing       test/timefp       test/testpoly
test/gendparam   test/testbb     test/testhilbert  test/testsig
test/checkfp
)
    add_executable(${bin} ${bin}.c)
    target_link_libraries(${bin} pbc ${GMP_LIBRARY} m)
endforeach(bin)
foreach(bin
    test/testbb test/testbls test/testibs
)
    target_link_libraries(${bin} pbc_sig)
endforeach(bin)

set(cppflags "-Iinclude")
set(warnflags "-Wall -W -Wfloat-equal -Wendif-labels -Wshadow  \\
    -Wpointer-arith -Wcast-align -Wstrict-prototypes \\
    -Wredundant-decls"
)
set(excessivewarnflags "-std=c99 -pedantic")
set(optflags "-O3 -ffast-math -fomit-frame-pointer \\
    -Wall -W -Wfloat-equal -Wendif-labels -Wshadow \\
    -Wpointer-arith -Wcast-align -Wstrict-prototypes  \\
    -Wredundant-decls"
)
set(CMAKE_C_FLAGS "-pipe ${cppflags} ${warnflags} ${optflags}")

#set(CMAKE_COMPILER_IS_GNUCXX 1)
#set(CMAKE_CXX_COMPILER "i586-mingw32msvc-c++")
#set(CMAKE_C_COMPILER "i586-mingw32msvc-gcc")
#set(CMAKE_AR "i586-mingw32msvc-ar")
#set(CMAKE_RANLIB "i586-mingw32msvc-ranlib")
#set(cppflags "${cppflags} -I/tmp/cross/include")
#set(CMAKE_C_FLAGS "-pipe ${cppflags} ${warnflags} ${optflags}")
