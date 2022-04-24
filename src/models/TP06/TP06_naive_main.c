#include <stdlib.h>

#include "TP06_naive.h"
#include "bench.h"


int main(int argc, char *argv[])
{
    struct options opt;
    bench_parse_args(argc, argv, &opt);
    bench(&model_TP06_naive, &opt);

    return EXIT_SUCCESS;
}
