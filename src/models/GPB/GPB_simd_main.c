#include <stdlib.h>

#include "../../bench.h"
#include "GPB_simd.h"


int main(int argc, char *argv[])
{
    struct options opt;
    bench_parse_args(argc, argv, &opt);
    bench(&model_GPB_simd, &opt);

    return EXIT_SUCCESS;
}
