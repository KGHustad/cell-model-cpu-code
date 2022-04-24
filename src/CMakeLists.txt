cmake_minimum_required (VERSION 3.9...3.17 FATAL_ERROR)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project("Cell model SIMD"
    VERSION 0.1
    LANGUAGES C CXX
)

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)

include(BuildType)
include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)
include(CheckIPOSupported)


option(BUILD_ACCURACY_CHECKS "Build programs to check accuracy. Requires the MPFR library." OFF)
option(HINT_CLANG_SIMD "Pass loop vectorisation pragma to clang" OFF)
option(HINT_OMP_SIMD "Pass simd directive to OpenMP" ON)
option(USE_FAST_MATH "Build with fast math enabled" ON)
option(USE_EXPM1 "Use expm1(x) where applicable instead of (exp(x) - 1). Increases accuracy, but may reduce performance." OFF)
option(USE_LTO "Build with link-time optimisation" ON)
option(WITH_OPENMP "Build with OpenMP" ON)

option(SET_MARCH "Pass machine-specific optimisation flags" ON)
option(TARGET_CASCADE_LAKE "Optimise for Intel Cascade Lake architecture" OFF)
option(TARGET_KNL "Optimise for Intel Knights Landing (KNL) architecture" OFF)
option(TARGET_SKYLAKE_256 "Optimise for Intel Skylake architecture using 256-bit instructions" OFF)
option(TARGET_SKYLAKE_512 "Optimise for Intel Skylake architecture using 512-bit instructions" OFF)

check_ipo_supported(RESULT ipo_is_supported)
if(USE_LTO AND ipo_is_supported)
    #message(STATUS "IPO is supported")
    # Enable inter-procedural optimisations (a.k.a. link-time optimisations)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    #set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)
endif()

if(HINT_CLANG_SIMD)
    add_definitions("-DHINT_CLANG_SIMD")
endif()
if(HINT_OMP_SIMD)
    add_definitions("-DHINT_OMP_SIMD")
endif()
if(LUT_INTERPOLATION_USE_NN)
    add_definitions("-DUSE_NN_LUT")
endif()
if(USE_EXPM1)
    add_definitions("-DUSE_EXPM1")
endif()

if(BUILD_ACCURACY_CHECKS)
    find_package(MPFR REQUIRED)
    if(NOT MPFR_FOUND)
        message(FATAL_ERROR "Could not find MPFR. Disabling BUILD_ACCURACY_CHECKS")
        set(BUILD_ACCURACY_CHECKS OFF CACHE BOOL "Build programs to check accuracy. Requires the MPFR library." FORCE)
    endif()
endif()

if(WITH_OPENMP)
    find_package(OpenMP)
    if(NOT OpenMP_FOUND)
        message(WARNING "OpenMP was not found. Disabling WITH_OPENMP")
        set(WITH_OPENMP OFF CACHE BOOL "Build with OpenMP" FORCE)
    endif()
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

message(STATUS "CPU architecture: ${CMAKE_SYSTEM_PROCESSOR}")


list(APPEND FAST_MATH_GNU_STYLE "GNU" "Clang" "AppleClang" "ArmClang" "FujitsuClang")
list(APPEND FAST_MATH_INTEL_STYLE "Intel")
set(FAST_MATH_FLAG_C "")
set(NO_FAST_MATH_FLAG_C "")
if(USE_FAST_MATH)
    if("${CMAKE_C_COMPILER_ID}" IN_LIST FAST_MATH_INTEL_STYLE)
        set(FAST_MATH_FLAG_C "-fp-model fast=2")
        #set(NO_FAST_MATH_FLAG_C "") # FIXME
    elseif("${CMAKE_C_COMPILER_ID}" IN_LIST FAST_MATH_GNU_STYLE)
        set(FAST_MATH_FLAG_C "-ffast-math")
        set(NO_FAST_MATH_FLAG_C "-fno-fast-math")
    else()
        message(WARNING "Fast math will not be enabled because the required flag for the detected compiler '${CMAKE_C_COMPILER_ID}' is unknown.")
    endif()
endif()

list(APPEND MARCH_GNU_STYLE "GNU" "Clang" "AppleClang" "Intel")
#list(APPEND MARCH_UNSUPPORTED "ARMClang")

set(MARCH_FLAG "")
message(STATUS "CMAKE_C_COMPILER_ID: ${CMAKE_C_COMPILER_ID}")
message(STATUS "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER}")

if(NOT SET_MARCH)
    # don't do anything
elseif(("${CMAKE_C_COMPILER_ID}" STREQUAL "ARMClang") OR (CMAKE_C_COMPILER MATCHES ".*armclang$"))
    # armclang does not support -march flag
    #message(WARNING "armclang does not support -march flag")
    if(USE_LTO)
        message(WARNING "-flto=thin breaks vectorisation. Please set USE_LTO to OFF.")
    endif()
    set(MARCH_FLAG "-fsimdmath -fno-math-errno -march=armv8.1-a -mcpu=native -flto")
    #set(MARCH_FLAG "-fsimdmath -fno-math-errno -Rpass-analysis=loop-vectorize -Ofast -march=armv8.1-a -mcpu=native -DHINT_CLANG_SIMD -flto")
elseif(("${CMAKE_C_COMPILER_ID}" STREQUAL "AppleClang") AND
    ("${CMAKE_SYSTEM_PROCESSOR}"
    STREQUAL "arm64"))
    message(WARNING "AppleClang on arm64 does not support -march flag")
elseif(("${CMAKE_C_COMPILER_ID}" STREQUAL "FujitsuClang") AND
    ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "aarch64"))
    # tune for Wisteria
    set(MARCH_FLAG "-Ofast -msve-vector-bits=512 -ffj-interleave-loop-insns=4")
elseif(("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang") AND
    ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64"))
    if(USE_LTO)
        message(WARNING "-flto=thin breaks vectorisation. Please set USE_LTO to OFF.")
    endif()
    #set(MARCH_FLAG "-march=native -fveclib=SVML -lsvml")
    set(MARCH_FLAG "-march=native -fveclib=libmvec")
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    if(TARGET_KNL)
        set(MARCH_FLAG "-axMIC-AVX512")
    elseif(TARGET_CASCADE_LAKE)
        set(MARCH_FLAG "-axCASCADELAKE")
    elseif(TARGET_SKYLAKE_512)
        set(MARCH_FLAG "-axCORE-AVX512")
    elseif(TARGET_SKYLAKE_256)
        set(MARCH_FLAG "-axSKYLAKE")
    else()
        set(MARCH_FLAG "-march=native")
    endif()
elseif("${CMAKE_C_COMPILER_ID}" IN_LIST MARCH_GNU_STYLE)
    if(TARGET_KNL)
        set(MARCH_FLAG "-march=knl -mprefer-vector-width=512")
    elseif(TARGET_CASCADE_LAKE)
        set(MARCH_FLAG "-march=native -mprefer-vector-width=512")
    elseif(TARGET_SKYLAKE_512)
        set(MARCH_FLAG "-march=native -mprefer-vector-width=512")
    elseif(TARGET_SKYLAKE_256)
        set(MARCH_FLAG "-march=native -mprefer-vector-width=256")
    else()
        set(MARCH_FLAG "-march=native")
    endif()
else()
    message(WARNING "Tuning for the native architecture will not be enabled because the required flag for the detected compiler '${CMAKE_C_COMPILER_ID}' is unknown.")
endif()

set(DEFAULT_VECTOR_LENGTH "")
if(("${CMAKE_C_COMPILER_ID}" STREQUAL "FujitsuClang") AND
    ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "aarch64"))
    set(DEFAULT_VECTOR_LENGTH "32")
endif()

if(VECTOR_LENGTH)
    set(VECTOR_LENGTH "${VECTOR_LENGTH}" CACHE STRING "Number of elements in vectors")
    set(MARCH_FLAG "${MARCH_FLAG} -DVECTOR_LENGTH=${VECTOR_LENGTH}")
else()
    set(VECTOR_LENGTH "${DEFAULT_VECTOR_LENGTH}" CACHE STRING "Number of elements in vectors" FORCE)
endif()

mark_as_advanced(CLEAR VECTOR_LENGTH)
if(VECTOR_LENGTH)
    message(STATUS "Using vector length ${VECTOR_LENGTH}")
endif()

set(CMAKE_C_FLAGS_DEBUG "-g -Wall -DDEBUG ${MARCH_FLAG}")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -DDEBUG ${MARCH_FLAG}")
check_c_compiler_flag("-Og" CC_SUPPORTS_OG)
check_cxx_compiler_flag("-Og" CXX_SUPPORTS_OG)
if(CC_SUPPORTS_OG)
    set(CMAKE_C_FLAGS_DEBUG "-Og ${CMAKE_C_FLAGS_DEBUG}")
endif()
if(CXX_SUPPORTS_OG)
    set(CMAKE_CXX_FLAGS_DEBUG "-Og ${CMAKE_CXX_FLAGS_DEBUG}")
endif()

set(CMAKE_C_FLAGS_MINSIZERELEASE "-Os ${MARCH_FLAG}")
set(CMAKE_CXX_FLAGS_MINSIZERELEASE "-Os ${MARCH_FLAG}")

set(CMAKE_C_FLAGS_RELEASE "-O3 ${MARCH_FLAG}")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 ${MARCH_FLAG}")

set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELEASE} -g")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -g")


include_directories(${PROJECT_SOURCE_DIR})

set(COMMON_BENCH_SOURCES bench.c bench_util.c common.c rng.cpp schemes.c)
set(COMMON_LUT_SOURCES
    lut/LUT.cpp
    lut/ExpressionLUT.cpp
    lut/LinearInterpolationLUT.cpp
)

# FAST_MATH_FLAG_C may contain multiple flags and must be passed as a list
string(REPLACE " " ";" FAST_MATH_FLAG_C_LIST "${FAST_MATH_FLAG_C}")

# define an interface to simplify linking
add_library(interface_bench INTERFACE)
target_compile_options(interface_bench INTERFACE ${FAST_MATH_FLAG_C_LIST})
target_link_options(interface_bench INTERFACE ${FAST_MATH_FLAG_C_LIST})
target_link_libraries(interface_bench INTERFACE m)
if(WITH_OPENMP AND OpenMP_C_FOUND)
    target_link_libraries(interface_bench INTERFACE OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
endif()

add_library(interface_bench_lut INTERFACE)
target_compile_options(interface_bench_lut INTERFACE ${FAST_MATH_FLAG_C_LIST})
target_link_options(interface_bench_lut INTERFACE ${FAST_MATH_FLAG_C_LIST})
if(WITH_OPENMP AND OpenMP_C_FOUND)
    target_link_libraries(interface_bench_lut INTERFACE OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
endif()


add_executable(math_bench_exp
    math_lib_bench/bench_exp.c
)
add_executable(math_bench_expm1
    math_lib_bench/bench_expm1.c
)
add_executable(math_bench_log
    math_lib_bench/bench_log.c
)
add_executable(math_bench_pow
    math_lib_bench/bench_pow.c
)
target_link_libraries(math_bench_exp PUBLIC interface_bench)
target_link_libraries(math_bench_expm1 PUBLIC interface_bench)
target_link_libraries(math_bench_log PUBLIC interface_bench)
target_link_libraries(math_bench_pow PUBLIC interface_bench)

if("${CMAKE_C_COMPILER_ID}" STREQUAL "AppleClang")
    find_library(ACCELERATE_LIBRARY Accelerate)
    target_link_libraries(math_bench_exp PUBLIC ${ACCELERATE_LIBRARY})
    target_link_libraries(math_bench_expm1 PUBLIC ${ACCELERATE_LIBRARY})
    target_link_libraries(math_bench_log PUBLIC ${ACCELERATE_LIBRARY})
    target_link_libraries(math_bench_pow PUBLIC ${ACCELERATE_LIBRARY})
endif()

add_executable(bench_TP06_naive
    ${COMMON_BENCH_SOURCES}
    models/TP06/TP06_naive.c
    models/TP06/TP06_naive_main.c
)
add_executable(bench_TP06_simd
    ${COMMON_BENCH_SOURCES}
    models/TP06/TP06_simd.c
    models/TP06/TP06_simd_main.c
)
add_executable(bench_TP06_nosimd
    ${COMMON_BENCH_SOURCES}
    models/TP06/TP06_nosimd.c
    models/TP06/TP06_nosimd_main.c
)
add_executable(compare_TP06_naive
    compare_naive_to_vanilla.c
    models/TP06/compare_TP06_naive.c
    models/TP06/TP06_naive.c
    models/TP06/TP06_simd.c
)
add_executable(compare_TP06_lut
    ${COMMON_LUT_SOURCES}
    compare_naive_to_lut.cpp
    models/TP06/compare_TP06_lut.cpp
    models/TP06/TP06_naive.c
    models/TP06/TP06_lut.cpp
)
add_executable(compare_TP06_simd_lut
    ${COMMON_LUT_SOURCES}
    compare_naive_to_lut.cpp
    models/TP06/compare_TP06_simd_lut.cpp
    models/TP06/TP06_naive.c
    models/TP06/TP06_simd_lut.cpp
)

add_executable(bench_JT21_naive
    ${COMMON_BENCH_SOURCES}
    models/JT21/JT21_naive.c
    models/JT21/JT21_naive_main.c
)
add_executable(bench_JT21_simd
    ${COMMON_BENCH_SOURCES}
    models/JT21/JT21_simd.c
    models/JT21/JT21_simd_main.c
)
add_executable(bench_JT21_nosimd
    ${COMMON_BENCH_SOURCES}
    models/JT21/JT21_nosimd.c
    models/JT21/JT21_nosimd_main.c
)
add_executable(compare_JT21_naive
    compare_naive_to_vanilla.c
    models/JT21/compare_JT21_naive.c
    models/JT21/JT21_naive.c
    models/JT21/JT21_simd.c
)
add_executable(bench_JT21_multistep_cell_time
    ${COMMON_BENCH_SOURCES}
    models/JT21/JT21_multistep_cell_time.c
    models/JT21/JT21_multistep_cell_time_main.c
)
add_executable(bench_JT21_multistep_cell_time_nosimd
    ${COMMON_BENCH_SOURCES}
    models/JT21/JT21_multistep_cell_time_nosimd.c
    models/JT21/JT21_multistep_cell_time_nosimd_main.c
)
add_executable(bench_JT21_multistep_cell_time_cell
    ${COMMON_BENCH_SOURCES}
    models/JT21/JT21_multistep_cell_time_cell.c
    models/JT21/JT21_multistep_cell_time_cell_main.c
)
add_executable(bench_JT21_multistep_cell_time_cell_nosimd
    ${COMMON_BENCH_SOURCES}
    models/JT21/JT21_multistep_cell_time_cell_nosimd.c
    models/JT21/JT21_multistep_cell_time_cell_nosimd_main.c
)
add_executable(bench_JT21_multistep_time_cell
    ${COMMON_BENCH_SOURCES}
    models/JT21/JT21_multistep_time_cell.c
    models/JT21/JT21_multistep_time_cell_main.c
)
add_executable(bench_JT21_multistep_time_cell_nosimd
    ${COMMON_BENCH_SOURCES}
    models/JT21/JT21_multistep_time_cell_nosimd.c
    models/JT21/JT21_multistep_time_cell_nosimd_main.c
)

add_executable(bench_GPB_naive
    ${COMMON_BENCH_SOURCES}
    models/GPB/GPB_naive.c
    models/GPB/GPB_naive_main.c
)
add_executable(bench_GPB_simd
    ${COMMON_BENCH_SOURCES}
    models/GPB/GPB_simd.c
    models/GPB/GPB_simd_main.c
)
add_executable(bench_GPB_nosimd
    ${COMMON_BENCH_SOURCES}
    models/GPB/GPB_nosimd.c
    models/GPB/GPB_nosimd_main.c
)
add_executable(compare_GPB_naive
    compare_naive_to_vanilla.c
    models/GPB/compare_GPB_naive.c
    models/GPB/GPB_naive.c
    models/GPB/GPB_simd.c
)

add_executable(TP06_dump_trace
    models/TP06/TP06_dump_trace.c
    models/TP06/TP06_naive.c
    models/TP06/TP06_simd.c
    common.c
    dump_trace_cmdline.c
    schemes.c
    solve.c
    trace.c
)
add_executable(TP06_dump_trace_lut
    ${COMMON_LUT_SOURCES}
    models/TP06/TP06_dump_trace_lut.cpp
    models/TP06/TP06_lut.cpp
    common.c
    schemes.c
    solve_lut.cpp
    trace.c
)

add_executable(bench_TP06_simd_lut
    ${COMMON_BENCH_SOURCES}
    ${COMMON_LUT_SOURCES}
    bench_lut.cpp
    models/TP06/TP06_simd_lut.cpp
    models/TP06/TP06_simd_lut_main.cpp
)
add_executable(bench_TP06_lut
    ${COMMON_BENCH_SOURCES}
    ${COMMON_LUT_SOURCES}
    bench_lut.cpp
    models/TP06/TP06_lut.cpp
    models/TP06/TP06_lut_main.cpp
)


target_link_libraries(TP06_dump_trace PUBLIC interface_bench)
target_link_libraries(TP06_dump_trace_lut PUBLIC interface_bench)

target_link_libraries(bench_JT21_naive PUBLIC interface_bench)
target_link_libraries(bench_JT21_simd PUBLIC interface_bench)
target_link_libraries(bench_JT21_nosimd PUBLIC interface_bench)
target_link_libraries(compare_JT21_naive PUBLIC interface_bench)

target_link_libraries(bench_JT21_multistep_cell_time PUBLIC interface_bench)
target_link_libraries(bench_JT21_multistep_cell_time_nosimd PUBLIC interface_bench)
target_link_libraries(bench_JT21_multistep_cell_time_cell PUBLIC interface_bench)
target_link_libraries(bench_JT21_multistep_cell_time_cell_nosimd PUBLIC interface_bench)
target_link_libraries(bench_JT21_multistep_time_cell PUBLIC interface_bench)
target_link_libraries(bench_JT21_multistep_time_cell_nosimd PUBLIC interface_bench)

target_link_libraries(bench_TP06_naive PUBLIC interface_bench)
target_link_libraries(bench_TP06_simd PUBLIC interface_bench)
target_link_libraries(bench_TP06_nosimd PUBLIC interface_bench)
target_link_libraries(compare_TP06_naive PUBLIC interface_bench)

target_link_libraries(bench_GPB_naive PUBLIC interface_bench)
target_link_libraries(bench_GPB_simd PUBLIC interface_bench)
target_link_libraries(bench_GPB_nosimd PUBLIC interface_bench)
target_link_libraries(compare_GPB_naive PUBLIC interface_bench)

target_link_libraries(bench_TP06_lut PUBLIC interface_bench_lut)
target_link_libraries(bench_TP06_simd_lut PUBLIC interface_bench_lut)


# Executables for verification of accuracy
add_executable(math_accuracy_compute
    math_lib_accuracy/compute.c
    math_lib_accuracy/common.c
)
target_link_libraries(math_accuracy_compute PUBLIC interface_bench)


if(BUILD_ACCURACY_CHECKS)
    add_executable(math_accuracy_check
        math_lib_accuracy/check.c
        math_lib_accuracy/common.c
    )

    target_include_directories(math_accuracy_check PUBLIC ${MPFR_INCLUDE_DIRS})
    target_link_libraries(math_accuracy_check PUBLIC ${MPFR_LIBRARIES})
    target_link_libraries(math_accuracy_check PUBLIC m)
endif()
