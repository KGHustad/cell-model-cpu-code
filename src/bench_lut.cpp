#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <vector>

#include "bench_lut.hpp"
#include "rng.h"
#include "schemes.h"

template <class LUT_type>
void bench_lut_generic_singlerep(const struct cellmodel_lut *model, struct options *opt,
                                 LUT_type &lut_V, LUT_type &lut_Ca)
{
    int num_threads = detect_num_threads();
    printf("Using %d threads\n", num_threads);

    printf("LUT class: %s\n", lut_V.class_desc);

    unsigned long num_cells = opt->num_cells;
    int num_steps = opt->num_timesteps;
    double dt = opt->dt;

    unsigned int num_states = model->num_states;
    unsigned int num_parameters = model->num_parameters;

    size_t alignment_bytes = CELLMODEL_STATES_ALIGNMENT_BYTES;
    unsigned int padded_num_cells = (unsigned int) ceil_to_multiple_uint64(
            num_cells, alignment_bytes / sizeof(cellmodel_float_t));

    size_t states_size = padded_num_cells * num_states * sizeof(cellmodel_float_t);
    size_t parameters_size = num_parameters * sizeof(cellmodel_float_t);

    cellmodel_float_t *states = (cellmodel_float_t *) aligned_alloc(alignment_bytes, states_size);
    cellmodel_float_t *parameters = (cellmodel_float_t *) malloc(parameters_size);

    model->init_states(states, num_cells, padded_num_cells);
    model->init_parameters(parameters);

    int is_stimulated_index = model->parameter_index("is_stimulated");
    if (is_stimulated_index == -1) {
        printf("Warning: Cannot apply stimulus when model has no is_stimulated parameter\n");
    }

    int V_index = -1;
    std::vector<const char *> V_names{"V_m", "V"};
    for (const char *V_name : V_names) {
        V_index = model->state_index(V_name);
        if (V_index != -1) {
            break;
        }
    }
    assert(V_index != -1);

    long V_offset = ((long) V_index) * num_cells;
    switch (opt->v_var) {
    case V_NOISE: {
        double noise_amplitude = 5;
        printf("V variation: Add random noise with amplitude %g\n", noise_amplitude);
        add_random_noise(states + V_offset, num_cells, opt->seed, noise_amplitude);
        break;
    }

    case V_LINEAR: {
        double amplitude = 5;
        printf("V variation: Linear displacement with amplitude %g\n", amplitude);
        for (unsigned long i = 0; i < num_cells; i++) {
            double ratio = i / (double) num_cells;
            states[V_offset + i] += 2 * amplitude * ratio - amplitude;
        }
        break;
    }

    case V_HOMOGENEOUS:
    default:
        printf("V variation: None (homogeneous)\n");
        break;
    }


    struct timespec timestamp_pre, timestamp_post;

    double t = 0;
    //double dt = 0.0001;
    clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_pre);
    step_lut_func step = NULL;
    switch (opt->scheme) {
    case SCHEME_FE:
        step = model->step_FE;
        if (step == NULL) {
            printf("ERROR: Model has no step function for scheme FE\n");
        }
        break;
    case SCHEME_RL:
        step = model->step_RL;
        if (step == NULL) {
            printf("ERROR: Model has no step function for scheme RL\n");
        }
        break;
    case SCHEME_GRL1:
        step = model->step_GRL1;
        if (step == NULL) {
            printf("ERROR: Model has no step function for scheme GRL1\n");
        }
        break;
    default:
        printf("ERROR: Unknown scheme with value %d\n", opt->scheme);
        exit(EXIT_FAILURE);
        break;
    }
    assert(step != NULL);


    for (int i = 0; i < num_steps; i++) {
        if (is_stimulated_index != -1) {
            if (t >= 1 && t < 2) {
                parameters[is_stimulated_index] = 1;
            } else {
                parameters[is_stimulated_index] = 0;
            }
        }
        step(states, t, dt, parameters, num_cells, padded_num_cells, lut_V, lut_Ca);
        t += dt;
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_post);

    double time_elapsed = (timestamp_post.tv_sec - timestamp_pre.tv_sec)
                          + 1E-9 * (timestamp_post.tv_nsec - timestamp_pre.tv_nsec);
    double cellsteps_per_second = num_cells * num_steps / time_elapsed;

    size_t tau_cellstep = 2 * num_states * sizeof(cellmodel_float_t);
    size_t tau = num_cells * num_steps * tau_cellstep;
    double mu_min = tau / time_elapsed;

    printf("Completed in %g seconds. %.4g cell steps per second. (%lu cells, %d steps.) Memory bandwidth >= %.4g GB/s\n",
           time_elapsed, cellsteps_per_second, num_cells, num_steps, mu_min * 1E-9);


    printf("V = %.15g at t = %g\n", states[padded_num_cells * V_index], t);

    free(states);
    free(parameters);
}

template <class LUT_type>
void bench_lut_generic_reps(const struct cellmodel_lut *model, struct options *opt, LUT_type &lut_V,
                            LUT_type &lut_Ca)
{
    int num_threads = detect_num_threads();
    printf("Using %d threads\n", num_threads);

    printf("LUT class: %s\n", lut_V.class_desc);

    unsigned long num_cells = opt->num_cells;
    int num_steps = opt->num_timesteps;
    double dt = opt->dt;
    int repetitions = opt->num_repetitions;

    unsigned int num_states = model->num_states;
    unsigned int num_parameters = model->num_parameters;

    size_t alignment_bytes = CELLMODEL_STATES_ALIGNMENT_BYTES;
    unsigned int padded_num_cells = (unsigned int) ceil_to_multiple_uint64(
            num_cells, alignment_bytes / sizeof(cellmodel_float_t));

    size_t states_size = padded_num_cells * num_states * sizeof(cellmodel_float_t);
    size_t parameters_size = num_parameters * sizeof(cellmodel_float_t);

    step_lut_func step = NULL;
    switch (opt->scheme) {
    case SCHEME_FE:
        step = model->step_FE;
        if (step == NULL) {
            printf("ERROR: Model has no step function for scheme FE\n");
        }
        break;
    case SCHEME_RL:
        step = model->step_RL;
        if (step == NULL) {
            printf("ERROR: Model has no step function for scheme RL\n");
        }
        break;
    case SCHEME_GRL1:
        step = model->step_GRL1;
        if (step == NULL) {
            printf("ERROR: Model has no step function for scheme GRL1\n");
        }
        break;
    default:
        printf("ERROR: Unknown scheme with value %d\n", opt->scheme);
        exit(EXIT_FAILURE);
        break;
    }
    assert(step != NULL);

    cellmodel_float_t *states = (cellmodel_float_t *) aligned_alloc(alignment_bytes, states_size);
    cellmodel_float_t *parameters = (cellmodel_float_t *) malloc(parameters_size);

    int stim_start_index = model->parameter_index("stim_start");

    int is_stimulated_index = model->parameter_index("is_stimulated");
    if (stim_start_index == -1 && is_stimulated_index == -1) {
        printf("Warning: Cannot apply stimulus when model has no 'stim_start' or 'is_stimulated' parameter\n");
    }

    int V_index = -1;
    std::vector<const char *> V_names{"V_m", "V"};
    for (const char *V_name : V_names) {
        V_index = model->state_index(V_name);
        if (V_index != -1) {
            break;
        }
    }
    assert(V_index != -1);

    double *time_elapsed_arr = (double *) malloc(sizeof(double) * repetitions);
    double *cellsteps_per_second_arr = (double *) malloc(sizeof(double) * repetitions);
    double *mu_min_arr = (double *) malloc(sizeof(double) * repetitions);

    double t = -1;
    for (int r = 0; r < repetitions; r++) {
        model->init_states(states, num_cells, padded_num_cells);
        model->init_parameters(parameters);

        if (stim_start_index != -1) {
            parameters[stim_start_index] = 1;
        }

        long V_offset = ((long) V_index) * num_cells;
        switch (opt->v_var) {
        case V_NOISE: {
            double noise_amplitude = 5;
            if (r == 0)
                printf("V variation: Add random noise with amplitude %g\n", noise_amplitude);
            add_random_noise(states + V_offset, num_cells, opt->seed, noise_amplitude);
            break;
        }

        case V_LINEAR: {
            double amplitude = 5;
            if (r == 0)
                printf("V variation: Linear displacement with amplitude %g\n", amplitude);
            for (unsigned long i = 0; i < num_cells; i++) {
                double ratio = i / (double) num_cells;
                states[V_offset + i] += 2 * amplitude * ratio - amplitude;
            }
            break;
        }

        case V_HOMOGENEOUS:
        default:
            if (r == 0)
                printf("V variation: None (homogeneous)\n");
            break;
        }

        struct timespec timestamp_pre, timestamp_post;
        t = 0;
        clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_pre);
        for (int i = 0; i < num_steps; i++) {
            if (is_stimulated_index != -1) {
                if (t >= 1 && t < 2) {
                    parameters[is_stimulated_index] = 1;
                } else {
                    parameters[is_stimulated_index] = 0;
                }
            }
            step(states, t, dt, parameters, num_cells, padded_num_cells, lut_V, lut_Ca);
            t += dt;
        }
        clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_post);

        double time_elapsed = (timestamp_post.tv_sec - timestamp_pre.tv_sec)
                              + 1E-9 * (timestamp_post.tv_nsec - timestamp_pre.tv_nsec);
        double cellsteps_per_second = num_cells * num_steps / time_elapsed;

        size_t tau_cellstep = 2 * num_states * sizeof(cellmodel_float_t);
        size_t tau = num_cells * num_steps * tau_cellstep;
        double mu_min = tau / time_elapsed;

        time_elapsed_arr[r] = time_elapsed;
        cellsteps_per_second_arr[r] = cellsteps_per_second;
        mu_min_arr[r] = mu_min;

        if (r == 0) {
            printf("%4s  %13s  %12s  %12s\n", "Rep", "Sol. time (s)", "Throughput", "Eff. mem. BW");
        }
        printf("%4d  %13g  %12.4e  %12.4e\n", r + 1, time_elapsed, cellsteps_per_second, mu_min);
    }

    qsort(time_elapsed_arr, repetitions, sizeof(double), &compare_double);
    qsort(cellsteps_per_second_arr, repetitions, sizeof(double), &compare_double_inv);
    qsort(mu_min_arr, repetitions, sizeof(double), &compare_double_inv);

    // best performance has lowest index
    int min_ind = repetitions - 1;
    int max_ind = 0;
    int best_ind = max_ind;
    int median_ind = repetitions / 2;
    printf("Best performance: %g seconds elapsed. %.4g cell steps per second. Memory bandwidth >= %.4g GB/s\n",
           time_elapsed_arr[best_ind], cellsteps_per_second_arr[best_ind],
           mu_min_arr[best_ind] * 1E-9);
#if 0
    printf("Median performance: %g seconds elapsed. %.4g cell steps per second. Memory bandwidth >= %.4g GB/s\n",
           time_elapsed_arr[median_ind], cellsteps_per_second_arr[median_ind],
           mu_min_arr[median_ind] * 1E-9);
#endif


    double min_relative_diff = relative_diff(mu_min_arr[min_ind], mu_min_arr[median_ind]);
    double max_relative_diff = relative_diff(mu_min_arr[max_ind], mu_min_arr[median_ind]);
    printf("Spread: [%.3g%%, %.3g%%] (worst and best performer's distance from median)\n",
           1E2 * min_relative_diff, 1E2 * max_relative_diff);

    size_t scheme_str_size = 20;
    char scheme_str[scheme_str_size];
    get_scheme_str(scheme_str, scheme_str_size, opt->scheme);
    printf("%lu cells, %d steps, %s scheme\n", num_cells, num_steps, scheme_str);


    printf("V = %.15g at t = %g\n", states[padded_num_cells * V_index], t);

    free(states);
    free(parameters);

    free(time_elapsed_arr);
    free(cellsteps_per_second_arr);
    free(mu_min_arr);
}

template <class LUT_type>
void bench_lut_generic(const struct cellmodel_lut *model, struct options *opt, LUT_type &lut_V,
                       LUT_type &lut_Ca)
{
    if (opt->num_repetitions > 1) {
        bench_lut_generic_reps<LUT_type>(model, opt, lut_V, lut_Ca);
    } else {
        bench_lut_generic_singlerep<LUT_type>(model, opt, lut_V, lut_Ca);
    }
}

void bench_lut(const struct cellmodel_lut *model, struct options *opt, default_LUT_type &lut_V,
               default_LUT_type &lut_Ca)
{
    bench_lut_generic<default_LUT_type>(model, opt, lut_V, lut_Ca);
}
