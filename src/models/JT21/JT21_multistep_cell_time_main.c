#include <stdlib.h>

#include "JT21_multistep_cell_time.h"
#include "bench.h"


int main(int argc, char *argv[])
{
    struct options opt;
    bench_parse_args(argc, argv, &opt);
    bench(&model_JT21_multistep_cell_time, &opt);

    return EXIT_SUCCESS;
}
