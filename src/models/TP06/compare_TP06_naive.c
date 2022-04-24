#include "../../compare_naive_to_vanilla.h"

#include "TP06_naive.h"
#include "TP06_simd.h"


int main()
{
    long num_cells = (long) 8;
    int num_steps = 40000;
    int show_rrms = 1;
    compare_naive_to_vanilla(&model_TP06_naive, &model_TP06_simd, num_cells, num_steps, SCHEME_FE, show_rrms);
    compare_naive_to_vanilla(&model_TP06_naive, &model_TP06_simd, num_cells, num_steps, SCHEME_RL, show_rrms);
    compare_naive_to_vanilla(&model_TP06_naive, &model_TP06_simd, num_cells, num_steps, SCHEME_GRL1, show_rrms);
}
