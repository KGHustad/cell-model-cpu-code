#include <stdlib.h>

#include "JT21_simd.h"
#include "bench.h"


int main(int argc, char *argv[])
{
    struct options opt;
    bench_parse_args(argc, argv, &opt);
    bench(&model_JT21_simd, &opt);

    return EXIT_SUCCESS;
}
