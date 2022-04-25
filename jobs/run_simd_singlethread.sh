BUILD_DIR=$(pwd)
if [ $# -gt 0 ]; then
    if [ -d $1 ] && [ -d "$1/bin" ]; then
        BUILD_DIR=$1
    else
        echo "Please provide the path to a build directory. Got $1"
    fi
fi

export OMP_NUM_THREADS=1

ARGS="--v_variation noise --seed 2021 -N 11688851 -t 2 --repetitions 5 -s FE"

implementations="naive nosimd simd"
models="TP06 JT21 GPB"

for model in $models; do
    for imp in $implementations; do
        echo
        echo "$model $imp"
        time "${BUILD_DIR}/bin/bench_${model}_${imp}" ${ARGS}
        echo
    done
    echo
done


date

echo "All done"
