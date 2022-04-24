#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cellmodel.h"
#include "compare.h"
#include "math_macros.h"

static uint64_t ceil_div_uint64(uint64_t a, uint64_t b)
{
    return (a + b - 1) / b;
}

static uint64_t ceil_to_multiple_uint64(uint64_t a, uint64_t b)
{
    return ceil_div_uint64(a, b) * b;
}

void compare_naive_to_vanilla(const struct cellmodel *model_naive,
                              const struct cellmodel *model_vanilla, long num_cells, int num_steps,
                              scheme_type scheme, int show_rrms)
{
    assert(model_naive->num_states == model_vanilla->num_states);
    assert(model_naive->num_parameters == model_vanilla->num_parameters);
    assert(model_naive->layout == LAYOUT_ARRAY_OF_STRUCTS);
    assert(model_vanilla->layout == LAYOUT_STRUCT_OF_ARRAYS);


    unsigned int num_states = model_naive->num_states;
    unsigned int num_parameters = model_naive->num_parameters;

    size_t alignment_bytes = CELLMODEL_STATES_ALIGNMENT_BYTES;
    unsigned int padded_num_cells = (unsigned int) ceil_to_multiple_uint64(
            num_cells, alignment_bytes / sizeof(cellmodel_float_t));

    size_t states_size = padded_num_cells * num_states * sizeof(cellmodel_float_t);
    size_t parameters_size = num_parameters * sizeof(cellmodel_float_t);


    cellmodel_float_t *states_naive = aligned_alloc(alignment_bytes, states_size);
    cellmodel_float_t *states_vanilla = aligned_alloc(alignment_bytes, states_size);
    //cellmodel_float_t *parameters = aligned_alloc(alignment_bytes, parameters_size);
    cellmodel_float_t *parameters_naive = malloc(parameters_size);
    cellmodel_float_t *parameters_vanilla = malloc(parameters_size);

    model_naive->init_states(states_naive, num_cells, padded_num_cells);
    model_vanilla->init_states(states_vanilla, num_cells, padded_num_cells);

    model_naive->init_parameters(parameters_naive);
    model_vanilla->init_parameters(parameters_vanilla);
    memcpy(parameters_vanilla, parameters_naive, parameters_size);


    cellmodel_float_t *states_error = aligned_alloc(alignment_bytes, states_size);

    double t = 0;
    double dt = 0.0005;
    int is_stimulated_index = model_naive->parameter_index("is_stimulated");

    int stim_start_index = model_naive->parameter_index("stim_start");
    if (stim_start_index != -1) {
        parameters_naive[stim_start_index] = 1;
        parameters_vanilla[stim_start_index] = 1;
    }

    cellmodel_float_t *states_error_Linf_array = malloc(sizeof(cellmodel_float_t) * num_steps);

    step_func model_naive_step = NULL;
    step_func model_vanilla_step = NULL;
    switch (scheme) {
    case SCHEME_FE:
        model_naive_step = model_naive->step_FE;
        if (model_naive_step == NULL) {
            printf("ERROR: Model A has no step function for scheme FE\n");
        }
        model_vanilla_step = model_vanilla->step_FE;
        if (model_vanilla_step == NULL) {
            printf("ERROR: Model B has no step function for scheme FE\n");
        }
        break;
    case SCHEME_RL:
        model_naive_step = model_naive->step_RL;
        if (model_naive_step == NULL) {
            printf("ERROR: Model A has no step function for scheme RL\n");
        }
        model_vanilla_step = model_vanilla->step_RL;
        if (model_vanilla_step == NULL) {
            printf("ERROR: Model B has no step function for scheme RL\n");
        }
        break;
    case SCHEME_GRL1:
        model_naive_step = model_naive->step_GRL1;
        if (model_naive_step == NULL) {
            printf("ERROR: Model A has no step function for scheme GRL1\n");
        }
        model_vanilla_step = model_vanilla->step_GRL1;
        if (model_vanilla_step == NULL) {
            printf("ERROR: Model B has no step function for scheme GRL1\n");
        }
        break;
    default:
        printf("ERROR: Unknown scheme with value %d\n", scheme);
        exit(EXIT_FAILURE);
        break;
    }
    assert(model_naive_step != NULL);
    assert(model_vanilla_step != NULL);

    // set up arrays to compute RRMS per state variable
    double numer_rrms[num_states];
    double denom_rrms[num_states];
    for (int s = 0; s < num_states; s++) {
        numer_rrms[s] = 0;
        denom_rrms[s] = 0;
    }

    for (int i = 0; i < num_steps; i++) {
        if (is_stimulated_index != -1) {
            if (t >= 1 && t < 2) {
                parameters_naive[is_stimulated_index] = 1;
                parameters_vanilla[is_stimulated_index] = 1;
            } else {
                parameters_naive[is_stimulated_index] = 0;
                parameters_vanilla[is_stimulated_index] = 0;
            }
        }

        model_naive_step(states_naive, t, dt, parameters_naive, num_cells, padded_num_cells);
        model_vanilla_step(states_vanilla, t, dt, parameters_vanilla, num_cells, padded_num_cells);

        //#pragma omp parallel for
        cellmodel_float_t states_error_Linf = 0;
        int max_s = -1;
        //int max_c = -1;
        for (int s = 0; s < num_states; s++) {
            {
                int c = 0;
                int ind_naive = c * num_states + s;
                int ind_vanilla = s * padded_num_cells + c;
                double approx = states_vanilla[ind_vanilla];
                double exact = states_naive[ind_naive];
                numer_rrms[s] += SQUARE(approx - exact);
                denom_rrms[s] += SQUARE(exact);
            }
            for (int c = 0; c < num_cells; c++) {
                int ind_naive = c * num_states + s;
                int ind_vanilla = s * padded_num_cells + c;
                states_error[ind_vanilla] = states_naive[ind_naive] - states_vanilla[ind_vanilla];
                if (fabs(states_error[ind_vanilla]) > states_error_Linf) {
                    states_error_Linf = fabs(states_error[ind_vanilla]);
                    max_s = s;
                    //max_c = c;
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

    double rrms_a[num_states];
    double rrms_max = 0;
    if (show_rrms)
        printf("RRMS for each state\n");
    for (int s = 0; s < num_states; s++) {
        rrms_a[s] = sqrt(numer_rrms[s] / denom_rrms[s]);
        const char *name = "";
        if (model_vanilla->state_names) {
            name = model_vanilla->state_names[s];
        }
        if (show_rrms)
            printf("State %3d (%-12s)  %.5e\n", s, name, rrms_a[s]);
        if (rrms_a[s] > rrms_max) {
            rrms_max = rrms_a[s];
        }
    }
    if (show_rrms)
        printf("Max RRMS: %g\n", rrms_max);

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

    printf("Max L-inf norm: %g at timestep %d (t=%g ms) with scheme %s\n",
           states_error_Linf_array_max, states_error_Linf_array_max_pos,
           states_error_Linf_array_max_pos * dt, scheme_str);

    free(states_error_Linf_array);
    free(states_error);

    free(states_naive);
    free(states_vanilla);
    free(parameters_naive);
    free(parameters_vanilla);
}
