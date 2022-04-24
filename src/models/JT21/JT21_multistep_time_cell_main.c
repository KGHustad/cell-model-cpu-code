#include <stdlib.h>

#include "JT21_multistep_time_cell.h"
#include "bench.h"


int main(int argc, char *argv[])
{
    struct options opt;
    bench_parse_args(argc, argv, &opt);
    bench(&model_JT21_multistep_time_cell, &opt);

    return EXIT_SUCCESS;
}
