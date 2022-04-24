#include "GPB_automatic.c"

#include "GPB_naive.h"

void init_states_omp(double *__restrict states, const long num_cells, const long padded_num_cells)
{
    #pragma omp parallel for
    for (int i = 0; i < num_cells; i++) {
        init_state_values(states + i * NUM_STATES);
    }
}

void GPB_step_FE(double *__restrict states, const double t, const double dt,
                 const double *__restrict parameters, const long num_cells, long padded_num_cells)
{
    #pragma omp parallel for
    for (int i = 0; i < num_cells; i++) {
        step_FE_singlecell(states + i * NUM_STATES, t, dt, parameters);
    }
}

void GPB_step_GRL1(double *__restrict states, const double t, const double dt,
                 const double *__restrict parameters, const long num_cells, long padded_num_cells)
{
    #pragma omp parallel for
    for (int i = 0; i < num_cells; i++) {
        step_GRL1_singlecell(states + i * NUM_STATES, t, dt, parameters);
    }
}

// clang-format off
const struct cellmodel model_GPB_naive = {
        .init_states = &init_states_omp,
        .init_parameters = &init_parameters_values,
        .state_index = &state_index,
        .parameter_index = &parameter_index,
        .step_FE = &GPB_step_FE,
        .step_RL = NULL,
        .step_GRL1 = &GPB_step_GRL1,
        .num_states = NUM_STATES,
        .num_parameters = NUM_PARAMS,
        .layout = LAYOUT_ARRAY_OF_STRUCTS,
        .colour_sets = NULL
};
// clang-format on
