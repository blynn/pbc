cmake_minimum_required(VERSION 2.4)
project(pbc)
include_directories(./include)

# cross compiling for Windows using Debian's mingw package
#   cmake -D crossmingw:BOOL=1
# I keep my cross-compiled GMP in /home/ben/cross/gmp/
if(crossmingw)
    include_directories(/home/ben/cross/gmp/include)
    set(GMP_LIBRARY "/home/ben/cross/gmp/lib/libgmp.a" )
    set(CMAKE_EXECUTABLE_SUFFIX ".exe")
    set(CMAKE_COMPILER_IS_GNUCXX 1)
    set(CMAKE_CXX_COMPILER "i586-mingw32msvc-c++")
    set(CMAKE_C_COMPILER "i586-mingw32msvc-gcc")
    set(CMAKE_AR "i586-mingw32msvc-ar")
    set(CMAKE_RANLIB "i586-mingw32msvc-ranlib")
    set(cppflags "${cppflags} -I/tmp/cross/include")
    set(CMAKE_C_FLAGS "-pipe ${cppflags} ${warnflags} ${optflags}")
    set(CMAKE_SKIP_RPATH ON)
    set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
    set(get_time_src "misc/get_time.win32.c")
else(crossmingw)
    find_library(GMP_LIBRARY NAMES gmp PATHS /usr/lib /usr/local/lib)
    set(get_time_src "misc/get_time.c")
endif(crossmingw)

set(leaksrc misc/leak.c)
set(unused
foo
)
add_library(pbc
arith/fieldmpz.c
arith/montfp.c
arith/naivefp.c
arith/tinyfp.c arith/fastfp.c arith/fasterfp.c
arith/field.c
arith/fp.c
arith/fieldquadratic.c
arith/indexcalculus.c
arith/poly.c arith/random.c
ecc/curve.c ecc/singular.c ecc/pairing.c ecc/param.c
ecc/a_param.c ecc/d_param.c ecc/e_param.c ecc/f_param.c
ecc/hilbert.c ecc/mnt.c ecc/mpc.c
misc/darray.c misc/symtab.c
misc/parse.c
misc/fops.c misc/tracker.c
misc/extend_printf.c
${get_time_src}
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
test/testexp
test/timersa
test/yuanli
test/testindexcalculus
)
    add_executable(${bin} ${bin}.c)
    target_link_libraries(${bin} pbc ${GMP_LIBRARY} m)
endforeach(bin)
foreach(bin
    test/testbb test/testbls test/testibs
)
    target_link_libraries(${bin} pbc_sig)
endforeach(bin)

#add_executable(test/pbc test/pbc.c test/pbc_getline.c)

add_executable(test/pbc test/pbc.c test/pbc_getline.readline.c)
target_link_libraries(test/pbc readline)

target_link_libraries(test/pbc pbc ${GMP_LIBRARY} m)

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

# my laptop (IBM Thinkpad X31)
#   cmake -D x31:BOOL=1 .
# doesn't help much
set(x31flags "-march=pentium-m -mmmx -msse -msse2 -mfpmath=sse -malign-double -maccumulate-outgoing-args")

set(CMAKE_C_FLAGS "-pipe ${cppflags} ${warnflags} ${optflags}")

if(x31)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${x31flags}")
endif(x31)
