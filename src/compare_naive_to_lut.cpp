#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cellmodel.h"
#include "cellmodel_lut.hpp"
#include "compare_naive_to_lut.hpp"

static uint64_t ceil_div_uint64(uint64_t a, uint64_t b)
{
    return (a + b - 1) / b;
}

static uint64_t ceil_to_multiple_uint64(uint64_t a, uint64_t b)
{
    return ceil_div_uint64(a, b) * b;
}

void compare_naive_to_lut(const struct cellmodel *model_naive,
                          const struct cellmodel_lut *model_lut, long num_cells, int num_steps,
                          double V_step, double Ca_step, scheme_type scheme)
{
    assert(model_naive->num_states == model_lut->num_states);
    assert(model_naive->num_parameters == model_lut->num_parameters);
    assert(model_naive->layout == LAYOUT_ARRAY_OF_STRUCTS);
    //assert(model_lut->layout == LAYOUT_STRUCT_OF_ARRAYS);


    unsigned int num_states = model_naive->num_states;
    unsigned int num_parameters = model_naive->num_parameters;

    size_t alignment_bytes = CELLMODEL_STATES_ALIGNMENT_BYTES;
    unsigned int padded_num_cells = (unsigned int) ceil_to_multiple_uint64(
            num_cells, alignment_bytes / sizeof(cellmodel_float_t));

    size_t states_size = padded_num_cells * num_states * sizeof(cellmodel_float_t);
    size_t parameters_size = num_parameters * sizeof(cellmodel_float_t);


    cellmodel_float_t *states_naive =
            (cellmodel_float_t *) aligned_alloc(alignment_bytes, states_size);
    cellmodel_float_t *states_lut =
            (cellmodel_float_t *) aligned_alloc(alignment_bytes, states_size);
    //cellmodel_float_t *parameters = aligned_alloc(alignment_bytes, parameters_size);
    cellmodel_float_t *parameters_naive = (cellmodel_float_t *) malloc(parameters_size);
    cellmodel_float_t *parameters_lut = (cellmodel_float_t *) malloc(parameters_size);

    model_naive->init_states(states_naive, num_cells, padded_num_cells);
    model_lut->init_states(states_lut, num_cells, padded_num_cells);

    model_naive->init_parameters(parameters_naive);
    model_lut->init_parameters(parameters_lut);
    memcpy(parameters_lut, parameters_naive, parameters_size);


    cellmodel_float_t *states_error =
            (cellmodel_float_t *) aligned_alloc(alignment_bytes, states_size);

    const int num_V_names = 2;
    int V_index = -1;
    const char V_names[][10] = {"V_m", "V"};
    for (int i = 0; i < num_V_names; i++) {
        V_index = model_lut->state_index(V_names[i]);
        if (V_index != -1) {
            break;
        }
    }

    double V_min = -100;
    double V_max = 100;
    //double V_step = 0.5;

    double Ca_min = 0.00001;
    double Ca_max = 10;
    //double Ca_step = 0.00001;

    double t = 0;
    double dt = 0.0005;
    int is_stimulated_index = model_naive->parameter_index("is_stimulated");

    int stim_start_index = model_naive->parameter_index("stim_start");
    if (stim_start_index != -1) {
        parameters_naive[stim_start_index] = 1;
        parameters_lut[stim_start_index] = 1;
    }

    cellmodel_float_t *states_error_Linf_array =
            (cellmodel_float_t *) malloc(sizeof(cellmodel_float_t) * num_steps);

    step_func model_naive_step = NULL;
    step_lut_func model_lut_step = NULL;
    switch (scheme) {
    case SCHEME_FE:
        model_naive_step = model_naive->step_FE;
        if (model_naive_step == NULL) {
            printf("ERROR: Model A has no step function for scheme FE\n");
        }
        model_lut_step = model_lut->step_FE;
        if (model_lut_step == NULL) {
            printf("ERROR: Model B has no step function for scheme FE\n");
        }
        break;
    case SCHEME_RL:
        model_naive_step = model_naive->step_RL;
        if (model_naive_step == NULL) {
            printf("ERROR: Model A has no step function for scheme RL\n");
        }
        model_lut_step = model_lut->step_RL;
        if (model_lut_step == NULL) {
            printf("ERROR: Model B has no step function for scheme RL\n");
        }
        break;
    case SCHEME_GRL1:
        model_naive_step = model_naive->step_GRL1;
        if (model_naive_step == NULL) {
            printf("ERROR: Model A has no step function for scheme GRL1\n");
        }
        model_lut_step = model_lut->step_GRL1;
        if (model_lut_step == NULL) {
            printf("ERROR: Model B has no step function for scheme GRL1\n");
        }
        break;
    default:
        printf("ERROR: Unknown scheme with value %d\n", scheme);
        exit(EXIT_FAILURE);
        break;
    }
    assert(model_naive_step != NULL);
    assert(model_lut_step != NULL);

    const std::vector<univariate_func> expressions_V = *model_lut->expressions_V;
    const std::vector<univariate_func> expressions_Ca = *model_lut->expressions_Ca;
    default_LUT_type lut_V(V_min, V_max, V_step, expressions_V, expressions_V.size(), dt,
                           parameters_lut);
    default_LUT_type lut_Ca(Ca_min, Ca_max, Ca_step, expressions_Ca, expressions_Ca.size(), dt,
                            parameters_lut);

    float_t states_error_rrms_numer[num_states];
    float_t states_error_rrms_denom[num_states];
    memset(states_error_rrms_numer, 0, sizeof(states_error_rrms_numer));
    memset(states_error_rrms_denom, 0, sizeof(states_error_rrms_denom));
    for (int i = 0; i < num_steps; i++) {
        if (is_stimulated_index != -1) {
            if (t >= 1 && t < 2) {
                parameters_naive[is_stimulated_index] = 1;
                parameters_lut[is_stimulated_index] = 1;
            } else {
                parameters_naive[is_stimulated_index] = 0;
                parameters_lut[is_stimulated_index] = 0;
            }
        }

        model_naive_step(states_naive, t, dt, parameters_naive, num_cells, padded_num_cells);
        model_lut_step(states_lut, t, dt, parameters_lut, num_cells, padded_num_cells, lut_V,
                       lut_Ca);

        //#pragma omp parallel for
        cellmodel_float_t states_error_Linf = 0;
        int max_s = -1;
        int max_c = -1;
        for (int s = 0; s < num_states; s++) {
            for (int c = 0; c < num_cells; c++) {
                int ind_naive = c * num_states + s;
                int ind_lut = s * padded_num_cells + c;
                float_t approx = states_lut[ind_lut];
                float_t exact = states_naive[ind_naive];
                float_t error = approx - exact;
                states_error[ind_lut] = error;
                if (fabs(states_error[ind_lut]) > states_error_Linf) {
                    states_error_Linf = fabs(states_error[ind_lut]);
                    max_s = s;
                    max_c = c;
                }

                if (c == 0) {
                    double exact = states_naive[ind_naive];
                    states_error_rrms_numer[s] += error * error;
                    states_error_rrms_denom[s] += exact * exact;
                }
            }
        }

        states_error_Linf_array[i] = states_error_Linf;

        if (num_steps <= 20) {
            printf("Step %2d. L-inf norm: %g at state %d\n", i, states_error_Linf, max_s);
        }

        t += dt;
    }

    cellmodel_float_t states_error_Linf_array_max = 0;
    int states_error_Linf_array_max_pos = -1;
    for (int i = 0; i < num_steps; i++) {
        if (states_error_Linf_array[i] > states_error_Linf_array_max) {
            states_error_Linf_array_max = states_error_Linf_array[i];
            states_error_Linf_array_max_pos = i;
        }
    }

    char scheme_str[100];
    switch (scheme) {
    case SCHEME_FE:
        sprintf(scheme_str, "FE");
        break;
    case SCHEME_RL:
        sprintf(scheme_str, "RL");
        break;
    case SCHEME_GRL1:
        sprintf(scheme_str, "GRL1");
        break;
    default:
        sprintf(scheme_str, "unknown");
    }

    double RRMS_V = sqrt(states_error_rrms_numer[V_index] / states_error_rrms_denom[V_index]);
    printf("RRMS for V: %-12g  ", RRMS_V);
    printf("Max L-inf norm: %g at timestep %d (t=%g ms) with scheme %s\n",
           states_error_Linf_array_max, states_error_Linf_array_max_pos,
           states_error_Linf_array_max_pos * dt, scheme_str);

    free(states_error_Linf_array);
    free(states_error);

    free(states_naive);
    free(states_lut);
    free(parameters_naive);
    free(parameters_lut);
}
