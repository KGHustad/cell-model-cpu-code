#include <stdio.h>
#include <stdlib.h>

#include <array>
#include <vector>

#include "TP06_simd_lut.hpp"
#include "bench_lut.hpp"
#include "lut/LUT.hpp"
#include "lut/LinearInterpolationLUT.hpp"


int main(int argc, char *argv[])
{
    const std::vector<univariate_func> expressions_V = *model_TP06_simd_lut.expressions_V;
    const std::vector<univariate_func> expressions_Ca = *model_TP06_simd_lut.expressions_Ca;

    printf("Num expressions (V) : %lu\n", expressions_V.size());
    printf("Num expressions (Ca): %lu\n", expressions_Ca.size());

    double V_min = -100;
    double V_max = 100;
    double V_step = 0.05;

    double Ca_min = 0.00001;
    double Ca_max = 10;
    double Ca_step = 0.00001;

    struct options opt;

    opt.V_step = V_step;
    opt.Ca_step = Ca_step;

    bench_parse_args(argc, argv, &opt);
    V_step = opt.V_step;
    Ca_step = opt.Ca_step;

    {
        double parameters[model_TP06_simd_lut.num_parameters];
        model_TP06_simd_lut.init_parameters(parameters);
        default_LUT_type lut_V(V_min, V_max, V_step, expressions_V, expressions_V.size(), opt.dt,
                               parameters);
        default_LUT_type lut_Ca(Ca_min, Ca_max, Ca_step, expressions_Ca, expressions_Ca.size(),
                                opt.dt, parameters);

        printf("V LUT memory footprint: %lu bytes\n", lut_V.get_size());
        printf("V LUT range [%g, %g] with step size %g\n", V_min, V_max, V_step);
        printf("Ca LUT memory footprint: %lu bytes\n", lut_Ca.get_size());
        printf("Ca LUT range [%g, %g] with step size %g\n", Ca_min, Ca_max, Ca_step);

        bench_lut(&model_TP06_simd_lut, &opt, lut_V, lut_Ca);
    }

    return EXIT_SUCCESS;
}
