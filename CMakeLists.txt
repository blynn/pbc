cmake_minimum_required(VERSION 2.4)
project(pbc)
set(CMAKE_MODULE_PATH .)
include_directories(. ./include)

CONFIGURE_FILE(config.h.cmake config.h)

set(examplebin
example/bls
example/bonehfranklin
example/hess
example/joux
example/paterson
example/yuanli
example/zhangkim
example/zss
)

set(goodbin
pbc/pbc
gen/gena1param
gen/genaparam
gen/gendparam
gen/geneparam
gen/genfparam
gen/gengparam
gen/hilbertpoly
gen/listmnt
gen/listfreeman
benchmark/benchmark
benchmark/timersa
benchmark/ellnet
)

set(buildbin ${goodbin} ${examplebin})

# cross compiling for Windows using Debian's mingw package
#   cmake -D crossmingw:BOOL=1
# I keep my cross-compiled GMP in ~/cross/gmp/
if(crossmingw)
    set(GMP_INCLUDE_DIR $ENV{HOME}/cross/gmp/include)
    set(GMP_LIBRARY $ENV{HOME}/cross/gmp/lib/libgmp.a)
    set(CMAKE_EXECUTABLE_SUFFIX ".exe")
    set(CMAKE_COMPILER_IS_GNUCXX 1)
    set(CMAKE_CXX_COMPILER "i586-mingw32msvc-c++")
    set(CMAKE_C_COMPILER "i586-mingw32msvc-gcc")
    set(CMAKE_AR "i586-mingw32msvc-ar")
    set(CMAKE_RANLIB "i586-mingw32msvc-ranlib")
    set(CMAKE_SKIP_RPATH ON)
    set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
    set(get_time_src "misc/get_time.win32.c")
else(crossmingw)
    find_package(GMP REQUIRED)
    find_package(Readline)
    set(get_time_src "misc/get_time.c")
endif(crossmingw)

if (READLINE_FOUND)
    add_executable(pbc/pbc pbc/pbc.c pbc/pbc_getline.readline.c)
    target_link_libraries(pbc/pbc readline)
else (READLINE_FOUND)
    add_executable(pbc/pbc pbc/pbc.c pbc/pbc_getline.c)
endif (READLINE_FOUND)

include_directories(${GMP_INCLUDE_DIR})

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
ecc/a_param.c ecc/d_param.c ecc/e_param.c ecc/f_param.c ecc/g_param.c
ecc/hilbert.c ecc/mnt.c ecc/mpc.c
misc/pbc_assert.c
misc/darray.c misc/symtab.c
misc/parse.c
misc/fops.c misc/tracker.c
misc/extend_printf.c
misc/memory.c
${get_time_src}
)

foreach(bin ${buildbin})
    add_executable(${bin} ${bin}.c)
    target_link_libraries(${bin} pbc ${GMP_LIBRARY} m)
endforeach(bin)

target_link_libraries(pbc/pbc pbc ${GMP_LIBRARY} m)

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

set(CMAKE_C_FLAGS "-pipe ${warnflags} ${optflags}")

if(x31)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${x31flags}")
endif(x31)

# e.g: cmake -D CMAKE_INSTALL_PREFIX=~/.local .
# make install

install(TARGETS pbc pbc LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib)

set(pbc_header
include/pbc_a1_param.h
include/pbc_a_param.h
include/pbc_assert.h
include/pbc_curve.h
include/pbc_darray.h
include/pbc_d_param.h
include/pbc_e_param.h
include/pbc_field.h
include/pbc_fieldmpz.h
include/pbc_fieldquadratic.h
include/pbc_fops.h
include/pbc_f_param.h
include/pbc_fp.h
include/pbc_g_param.h
include/pbc.h
include/pbc_hilbert.h
include/pbc_memory.h
include/pbc_mnt.h
include/pbc_mpc.h
include/pbc_pairing.h
include/pbc_param.h
include/pbc_parse.h
include/pbc_poly.h
include/pbc_random.h
include/pbc_singular.h
include/pbc_symtab.h
include/pbc_time.h
include/pbc_tracker.h
include/pbc_utils.h
)

install(FILES ${pbc_header} DESTINATION include/pbc)
