#include <stdlib.h>

#include "bench.h"
#include "JT21_naive.h"


int main(int argc, char *argv[])
{
    struct options opt;
    bench_parse_args(argc, argv, &opt);
    bench(&model_JT21_naive, &opt);

    return EXIT_SUCCESS;
}
