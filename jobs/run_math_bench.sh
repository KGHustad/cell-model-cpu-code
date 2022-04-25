BUILD_DIR=$(pwd)
if [ $# -gt 0 ]; then
    if [ -d $1 ] && [ -d "$1/bin" ]; then
        BUILD_DIR=$1
    else
        echo "Please provide the path to a build directory. Got $1"
    fi
fi

export OMP_NUM_THREADS=1

math_functions="exp expm1 log pow"
N=30000
repetitions=20000

for math_func in ${math_functions}; do
    echo
    echo "${math_func}"
    time "${BUILD_DIR}/bin/math_bench_${math_func}" ${N} ${repetitions}
    echo
done


date

echo "All done"
