#include "../../compare_naive_to_vanilla.h"

#include "GPB_naive.h"
#include "GPB_simd.h"


int main()
{
    long num_cells = (long) 8;
    int num_steps = 4000;
    int show_rrms = 0;
    compare_naive_to_vanilla(&model_GPB_naive, &model_GPB_simd, num_cells, num_steps, SCHEME_FE, show_rrms);
    compare_naive_to_vanilla(&model_GPB_naive, &model_GPB_simd, num_cells, num_steps, SCHEME_GRL1, show_rrms);
}
