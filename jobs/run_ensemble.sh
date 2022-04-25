BUILD_DIR=$(pwd)
if [ $# -gt 0 ]; then
    if [ -d $1 ] && [ -d "$1/bin" ]; then
        BUILD_DIR=$1
    else
        echo "Please provide the path to a build directory. Got $1"
        exit 1
    fi
fi

ARGS="--v_variation noise --seed 2021 -N 11688851 -t 100 --repetitions 5 --multistep"

implementations="cell_time cell_time_nosimd cell_time_cell cell_time_cell_nosimd time_cell time_cell_nosimd"

for imp in $implementations; do
    echo
    echo "$imp"
    echo
    time "${BUILD_DIR}/bin/bench_JT21_multistep_${imp}" ${ARGS} -s FE
    echo
    time "${BUILD_DIR}/bin/bench_JT21_multistep_${imp}" ${ARGS} -s GRL1
    echo
done
echo

date

echo "All done"
