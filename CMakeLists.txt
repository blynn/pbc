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

    add_executable(test/pbc test/pbc.c test/pbc_getline.c)

else(crossmingw)
    find_library(GMP_LIBRARY NAMES gmp PATHS /usr/lib /usr/local/lib)
    set(get_time_src "misc/get_time.c")

    add_executable(test/pbc test/pbc.c test/pbc_getline.readline.c)
    target_link_libraries(test/pbc readline)

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
misc/pbc_assert.c
misc/darray.c misc/symtab.c
misc/parse.c
misc/fops.c misc/tracker.c
misc/extend_printf.c
${get_time_src}
)
foreach(bin
test/19          test/geneparam  test/testibe
test/benchmark   test/genfparam  test/testsig
test/gena1param  test/listmnt    test/testfi       test/testpairing
test/genaparam   test/sing       test/timefp       test/testpoly
test/gendparam   test/testhilbert
test/checkfp
test/testexp
test/timersa
test/yuanli
test/testindexcalculus
)
    add_executable(${bin} ${bin}.c)
    target_link_libraries(${bin} pbc ${GMP_LIBRARY} m)
endforeach(bin)

target_link_libraries(test/pbc pbc ${GMP_LIBRARY} m)

set(cppflags "-Iinclude")
set(warnflags "-Wall -W -Wfloat-equal -Wendif-labels -Wshadow  \\
    -Wpointer-arith -Wcast-align -Wstrict-prototypes \\
    -Wredundant-decls"
)
set(excessivewarnflags "-std=c99 -pedantic")

if(debug)
set(optflags "-O0 -g")
else(debug)
set(optflags "-O3 -ffast-math -fomit-frame-pointer")
endif(debug)

# my laptop (IBM Thinkpad X31)
#   cmake -D x31:BOOL=1 .
# doesn't help much
set(x31flags "-march=pentium-m -mmmx -msse -msse2 -mfpmath=sse -malign-double -maccumulate-outgoing-args")

set(CMAKE_C_FLAGS "-pipe ${cppflags} ${warnflags} ${optflags}")

if(x31)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${x31flags}")
endif(x31)

# e.g: cmake -D CMAKE_INSTALL_PREFIX=~/.local .
# make install

install(TARGETS pbc pbc LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib)

set(pbc_header
include/a1_param.h
include/a_param.h
include/curve.h
include/darray.h
include/d_param.h
include/e_param.h
include/field.h
include/fieldmpz.h
include/fieldquadratic.h
include/fops.h
include/f_param.h
include/fp.h
include/get_time.h
include/hash.h
include/hilbert.h
include/mnt.h
include/mpc.h
include/pairing.h
include/param.h
include/parse.h
include/pbc_assert.h
include/pbc.h
include/poly.h
include/random.h
include/singular.h
include/symtab.h
include/tracker.h
include/utils.h
)

install(FILES ${pbc_header} DESTINATION include/pbc)
