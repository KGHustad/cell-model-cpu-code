#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cellmodel.h"
#include "compare.h"

static uint64_t ceil_div_uint64(uint64_t a, uint64_t b)
{
    return (a + b - 1) / b;
}

static uint64_t ceil_to_multiple_uint64(uint64_t a, uint64_t b)
{
    return ceil_div_uint64(a, b) * b;
}

void compare(const struct cellmodel *model_a, const struct cellmodel *model_b, long num_cells,
             int num_steps, scheme_type scheme)
{
    assert(model_a->num_states == model_b->num_states);
    assert(model_a->num_parameters == model_b->num_parameters);

    unsigned int num_states = model_a->num_states;
    unsigned int num_parameters = model_a->num_parameters;

    size_t alignment_bytes = CELLMODEL_STATES_ALIGNMENT_BYTES;
    unsigned int padded_num_cells = (unsigned int) ceil_to_multiple_uint64(
            num_cells, alignment_bytes / sizeof(cellmodel_float_t));

    size_t states_size = padded_num_cells * num_states * sizeof(cellmodel_float_t);
    //size_t parameters_size = padded_num_cells * num_parameters * sizeof(cellmodel_float_t);
    size_t parameters_size = num_parameters * sizeof(cellmodel_float_t);


    cellmodel_float_t *states_a = aligned_alloc(alignment_bytes, states_size);
    cellmodel_float_t *states_b = aligned_alloc(alignment_bytes, states_size);
    //cellmodel_float_t *parameters = aligned_alloc(alignment_bytes, parameters_size);
    cellmodel_float_t *parameters_a = malloc(parameters_size);
    cellmodel_float_t *parameters_b = malloc(parameters_size);

    model_a->init_states(states_a, num_cells, padded_num_cells);
    model_b->init_states(states_b, num_cells, padded_num_cells);
    memcpy(states_b, states_a, states_size);

    model_a->init_parameters(parameters_a);
    model_b->init_parameters(parameters_b);
    memcpy(parameters_b, parameters_a, parameters_size);


    cellmodel_float_t *states_error = aligned_alloc(alignment_bytes, states_size);

    double t = 0;
    double dt = 0.0005;
    int is_stimulated_index = model_a->parameter_index("is_stimulated");

    int stim_start_index = model_a->parameter_index("stim_start");
    if (stim_start_index != -1) {
        parameters_a[stim_start_index] = 1;
        parameters_b[stim_start_index] = 1;
    }

    cellmodel_float_t *states_error_Linf_array = malloc(sizeof(cellmodel_float_t) * num_steps);

    step_func model_a_step = NULL;
    step_func model_b_step = NULL;
    switch (scheme) {
    case SCHEME_FE:
        model_a_step = model_a->step_FE;
        if (model_a_step == NULL) {
            printf("ERROR: Model A has no step function for scheme FE\n");
        }
        model_b_step = model_b->step_FE;
        if (model_b_step == NULL) {
            printf("ERROR: Model B has no step function for scheme FE\n");
        }
        break;
    case SCHEME_RL:
        model_a_step = model_a->step_RL;
        if (model_a_step == NULL) {
            printf("ERROR: Model A has no step function for scheme RL\n");
        }
        model_b_step = model_b->step_RL;
        if (model_b_step == NULL) {
            printf("ERROR: Model B has no step function for scheme RL\n");
        }
        break;
    case SCHEME_GRL1:
        model_a_step = model_a->step_GRL1;
        if (model_a_step == NULL) {
            printf("ERROR: Model A has no step function for scheme GRL1\n");
        }
        model_b_step = model_b->step_GRL1;
        if (model_b_step == NULL) {
            printf("ERROR: Model B has no step function for scheme GRL1\n");
        }
        break;
    default:
        printf("ERROR: Unknown scheme with value %d\n", scheme);
        exit(EXIT_FAILURE);
        break;
    }
    assert(model_a_step != NULL);
    assert(model_b_step != NULL);


    for (int i = 0; i < num_steps; i++) {
        if (is_stimulated_index != -1) {
            if (t >= 1 && t < 2) {
                parameters_a[is_stimulated_index] = 1;
                parameters_b[is_stimulated_index] = 1;
            } else {
                parameters_a[is_stimulated_index] = 0;
                parameters_b[is_stimulated_index] = 0;
            }
        }

        model_a_step(states_a, t, dt, parameters_a, num_cells, padded_num_cells);
        model_b_step(states_b, t, dt, parameters_b, num_cells, padded_num_cells);

        //#pragma omp parallel for
        cellmodel_float_t states_error_Linf = 0;
        int max_s = -1;
        int max_c = -1;
        for (int s = 0; s < num_states; s++) {
            for (int c = 0; c < num_cells; c++) {
                int ind = s * padded_num_cells + c;
                states_error[ind] = states_a[ind] - states_b[ind];
                if (fabs(states_error[ind]) > states_error_Linf) {
                    states_error_Linf = fabs(states_error[ind]);
                    max_s = s;
                    max_c = c;
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

    printf("Max L-inf norm: %g at timestep %d with scheme %s\n", states_error_Linf_array_max,
           states_error_Linf_array_max_pos, scheme_str);

    free(states_error_Linf_array);
    free(states_error);

    free(states_a);
    free(states_b);
    free(parameters_a);
    free(parameters_b);
}
