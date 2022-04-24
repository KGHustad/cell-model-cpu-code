#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cellmodel_lut.hpp"
#include "schemes.h"
#include "solve_lut.hpp"

static uint64_t ceil_div_uint64(uint64_t a, uint64_t b)
{
    return (a + b - 1) / b;
}

static uint64_t ceil_to_multiple_uint64(uint64_t a, uint64_t b)
{
    return ceil_div_uint64(a, b) * b;
}

void solve_single_lut(const struct cellmodel_lut *model, double dt, double T_end, int store_period,
                      scheme_type scheme, cellmodel_float_t **V_trace_ptr,
                      cellmodel_float_t **t_trace_ptr, int *trace_length, default_LUT_type &lut_V,
                      default_LUT_type &lut_Ca)
{
    unsigned int num_states = model->num_states;
    unsigned int num_parameters = model->num_parameters;

    int num_cells = 1;
    size_t alignment_bytes = CELLMODEL_STATES_ALIGNMENT_BYTES;
    unsigned int padded_num_cells = (unsigned int) ceil_to_multiple_uint64(
            num_cells, alignment_bytes / sizeof(cellmodel_float_t));

    size_t states_size = padded_num_cells * num_states * sizeof(cellmodel_float_t);
    size_t parameters_size = num_parameters * sizeof(cellmodel_float_t);

    cellmodel_float_t *states = (cellmodel_float_t *) aligned_alloc(alignment_bytes, states_size);
    cellmodel_float_t *states_tmp1 = NULL;
    cellmodel_float_t *states_tmp2 = NULL;
    cellmodel_float_t *parameters = (cellmodel_float_t *) malloc(parameters_size);

    model->init_states(states, num_cells, padded_num_cells);
    model->init_parameters(parameters);

    long num_timesteps = lround(T_end / dt);
    T_end = num_timesteps * dt;

    int num_stored_timesteps = num_timesteps / store_period;
    num_stored_timesteps++;

    cellmodel_float_t *t_trace =
            (cellmodel_float_t *) malloc(num_stored_timesteps * sizeof(cellmodel_float_t));
    cellmodel_float_t *V_trace =
            (cellmodel_float_t *) malloc(num_stored_timesteps * sizeof(cellmodel_float_t));

    const int num_V_names = 2;
    int V_index = -1;
    const char V_names[][10] = {"V_m", "V"};
    for (int i = 0; i < num_V_names; i++) {
        V_index = model->state_index(V_names[i]);
        if (V_index != -1) {
            break;
        }
    }

    step_lut_func step = NULL;
    switch (scheme) {
    case SCHEME_FE:
        step = model->step_FE;
        if (step == NULL) {
            printf("ERROR: Model has no step function for scheme FE\n");
            exit(EXIT_FAILURE);
        }
        break;
    case SCHEME_RL:
        step = model->step_RL;
        if (step == NULL) {
            printf("ERROR: Model has no step function for scheme RL\n");
            exit(EXIT_FAILURE);
        }
        break;
    case SCHEME_GRL1:
        step = model->step_GRL1;
        if (step == NULL) {
            printf("ERROR: Model has no step function for scheme GRL1\n");
            exit(EXIT_FAILURE);
        }
        break;
    default:
        printf("ERROR: Unknown scheme with value %d\n", scheme);
        exit(EXIT_FAILURE);
        break;
    }
    //assert(step != NULL);

    // set stimulus parameters
    double eps = 1E-8;
    int is_stimulated_index = model->parameter_index("is_stimulated");

    int stim_start_index = model->parameter_index("stim_start");
    if (stim_start_index != -1) {
        parameters[stim_start_index] = 0;
    }

    //printf("num_timesteps: %ld\n", num_timesteps);
    double t = 0;
    int trace_pos = 0;
    t_trace[trace_pos] = t;
    unsigned int V_ind_2d = V_index * padded_num_cells;
    V_trace[trace_pos] = states[V_ind_2d];
    trace_pos++;
    for (long i = 1; i <= num_timesteps; i++) {
        if (is_stimulated_index != -1) {
            if (t >= 1 - eps && t < 2 - eps) {
                parameters[is_stimulated_index] = 1;
            } else {
                parameters[is_stimulated_index] = 0;
            }
        }
        step(states, t, dt, parameters, num_cells, padded_num_cells, lut_V, lut_Ca);
        t += dt;

        //printf("At t=%g\n", t);

        if ((i % store_period) == 0) {
            assert(trace_pos < num_stored_timesteps);
            t_trace[trace_pos] = t;
            V_trace[trace_pos] = states[V_ind_2d];
            trace_pos++;
        }
    }

    free(states);
    free(parameters);
    if (states_tmp1) {
        free(states_tmp1);
    }
    if (states_tmp2) {
        free(states_tmp2);
    }

    *V_trace_ptr = V_trace;
    *t_trace_ptr = t_trace;
    *trace_length = num_stored_timesteps;
}
