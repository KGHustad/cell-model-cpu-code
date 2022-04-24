
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "TP06_lut.hpp"
#include "schemes.h"
#include "solve_lut.hpp"
#include "trace.h"


int main()
{
    const std::vector<univariate_func> expressions_V = *model_TP06_lut.expressions_V;
    const std::vector<univariate_func> expressions_Ca = *model_TP06_lut.expressions_Ca;

    const struct cellmodel_lut *model = &model_TP06_lut;
    double T_end = 1000;
    double solve_dt = 1E-3;
    int store_period = 1;

    double V_min = -100;
    double V_max = 100;
    double V_step = 0.05;

    double Ca_min = 0.00001;
    double Ca_max = 10;
    double Ca_step = 0.00001;

#if 0
    T_end = 10;
    solve_dt = 0.1;
    store_period = 1;
#endif

    scheme_type scheme = SCHEME_RL;
    cellmodel_float_t *V_trace;
    cellmodel_float_t *t_trace;
    int trace_length;

    double parameters[model_TP06_lut.num_parameters];
    model_TP06_lut.init_parameters(parameters);
    default_LUT_type lut_V(V_min, V_max, V_step, expressions_V, expressions_V.size(), solve_dt,
                            parameters);
    default_LUT_type lut_Ca(Ca_min, Ca_max, Ca_step, expressions_Ca, expressions_Ca.size(),
                            solve_dt, parameters);

    long num_timesteps = lround(T_end / solve_dt);

    printf("Solving with %ld steps\n", num_timesteps);

    solve_single_lut(model, solve_dt, T_end, store_period, scheme, &V_trace, &t_trace, &trace_length, lut_V, lut_Ca);

    printf("Storing trace of length %d\n", trace_length);

    size_t scheme_str_size = 20;
    char scheme_str[scheme_str_size];
    get_scheme_str(scheme_str, scheme_str_size, scheme);

    char basename[200];
    sprintf(basename, "dat/lut/TP06_%s_%g", scheme_str, solve_dt);

    dump_trace(t_trace, V_trace, trace_length, solve_dt, scheme, basename);
}
