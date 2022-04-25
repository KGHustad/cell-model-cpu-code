# Cell model implementations optimised for execution on modern CPUs

Source code for the implementations discussed in the paper "Resource-efficient use of modern processor architectures for numerically solving cardiac ionic cell models" by Kristian Gregorius Hustad and Xing Cai.

## Dependencies
- [CMake](https://cmake.org/)
- C and C++ compilers (see compiler-specific details below)

## Building
```bash
build_dir="build"
rm -rf ${build_dir}
mkdir -p ${build_dir}
pushd ${build_dir}

cmake ../src
make -j8 # replace 8 by the number of parallel build jobs you want to run

popd
```

### Fujitsu compiler
CMake 3.21 or newer is required for proper support of the Fujitsu compiler in Clang mode.

### Controlling which compiler is used with CMake
- ARMclang: `export CC=armclang CXX=armclang++`
- Clang: `export CC=clang CXX=clang++`
- Fujitsu (Clang mode): `export CC="fcc -Nclang" CXX="FCC -Nclang"`
- GCC: `export CC=gcc CXX=g++`
- Intel: `export CC=icc CXX=icpc`

### SIMD hints
We use ifdefs to determine how SIMD hints should be provided.

For Clang-based compilers (i.e. bog-standard clang, ArmClang, Fujitsu compiler in Clang mode), `HINT_CLANG_SIMD` should be enabled.

For other compilers (e.g. GCC, Intel), `HINT_OMP_SIMD` should be enabled. (Note that `HINT_CLANG_SIMD` is checked before `HINT_OMP_SIMD`, so there is no harm in leaving `HINT_OMP_SIMD` enabled when using a Clang-based compiler.)
When using `HINT_OMP_SIMD`, `-DVECTOR_LENGTH=8` can be passed to CMake to extend the hint with an additional clause `simdlen(VECTOR_LENGTH)`. This is needed to ensure AVX-512 is used on Cascade Lake and Skylake (but not KNL) with the Intel compiler.

### Target-specific options
Only one of the following options should be set.
- `TARGET_CASCADE_LAKE`: For the Cascade Lake CPUs (found on Oakbridge)
- `TARGET_KNL`: For the Knights Landing CPU (found on Oakforest)
- `TARGET_SKYLAKE_256`: For the Skylake client CPUs (without AVX-512)
- `TARGET_SKYLAKE_512`: For the Skylake (and later *lake) server CPUs (with AVX-512)

If none of these options are set, the compiler is set to target the native CPU architecture (if supported).

__Note: The compiler used will affect which optimisation flags are passed by CMake.__

### Other options
- `USE_EXPM1`: Use the `expm1()` function to compute expressions on the form `(exp(x) - 1)`.
- `USE_FAST_MATH`: Pass flags to the compiler that increase math performance, potentially at the cost of a slight loss of accuracy.
- `USE_LTO`: Use link-time optimisation (LTO) / interprocedural optimisation (IPO)
- `VECTOR_LENGTH`: Override the vector length. (It is safe to set this to a multiple of the hardware vector length.)

## Running
We have provided scripts to run the same benchmarks that were used to generate the data listed in the following tables:
- Table 3: `jobs/run_math_bench.sh`
- Table 6: `jobs/run_simd_singlethread.sh`
- Table 7: `jobs/run_simd_multithread.sh`
- Table 8: `jobs/run_ensemble.sh`
- Table 10: `jobs/run_TP06_lut.sh`

These scripts take the build directory as an argument. Example usage:
```
bash jobs/run_math_TP06_lut.sh build
```
