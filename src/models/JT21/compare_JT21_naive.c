#include "../../compare_naive_to_vanilla.h"

#include "JT21_naive.h"
#include "JT21_simd.h"


int main()
{
    long num_cells = (long) 8;
    int num_steps = 4000;
    int show_rrms = 0;
    compare_naive_to_vanilla(&model_JT21_naive, &model_JT21_simd, num_cells, num_steps, SCHEME_FE, show_rrms);
    compare_naive_to_vanilla(&model_JT21_naive, &model_JT21_simd, num_cells, num_steps, SCHEME_GRL1, show_rrms);
}
