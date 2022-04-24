#include "TP06_automatic.c"

#include "TP06_naive.h"

void init_states_omp(double *__restrict states, const long num_cells, const long padded_num_cells)
{
    #pragma omp parallel for
    for (int i = 0; i < num_cells; i++) {
        init_state_values(states + i * NUM_STATES);
    }
}

void TP06_step_FE(double *__restrict states, const double t, const double dt,
                  const double *__restrict parameters, const long num_cells, long padded_num_cells)
{
    #pragma omp parallel for
    for (int i = 0; i < num_cells; i++) {
        step_FE_singlecell(states + i * NUM_STATES, t, dt, parameters);
    }
}

void TP06_step_GRL1(double *__restrict states, const double t, const double dt,
                    const double *__restrict parameters, const long num_cells,
                    long padded_num_cells)
{
    #pragma omp parallel for
    for (int i = 0; i < num_cells; i++) {
        step_GRL1_singlecell(states + i * NUM_STATES, t, dt, parameters);
    }
}

void TP06_step_RL(double *__restrict states, const double t, const double dt,
                  const double *__restrict parameters, const long num_cells, long padded_num_cells)
{
    #pragma omp parallel for
    for (int i = 0; i < num_cells; i++) {
        step_RL_singlecell(states + i * NUM_STATES, t, dt, parameters);
    }
}

// clang-format off
const struct cellmodel model_TP06_naive = {
        .init_states = &init_states_omp,
        .init_parameters = &init_parameters_values,
        .state_index = &state_index,
        .parameter_index = &parameter_index,
        .step_FE = &TP06_step_FE,
        .step_RL = &TP06_step_RL,
        .step_GRL1 = &TP06_step_GRL1,
        .num_states = NUM_STATES,
        .num_parameters = NUM_PARAMS,
        .layout = LAYOUT_ARRAY_OF_STRUCTS,
        .colour_sets = NULL
};
// clang-format on
