BUILD_DIR=$(pwd)
if [ $# -gt 0 ]; then
    if [ -d $1 ] && [ -d "$1/bin" ]; then
        BUILD_DIR=$1
    else
        echo "Please provide the path to a build directory. Got $1"
    fi
fi

ARGS="--v_variation noise --seed 2021 -N 11688851 -t 100 --repetitions 5 -s RL"

implementations="naive nosimd simd lut simd_lut"
model="TP06"

for imp in $implementations; do
    echo
    echo "$model $imp"
    time "${BUILD_DIR}/bin/bench_${model}_${imp}" ${ARGS}
    echo
done
echo


date

echo "All done"
