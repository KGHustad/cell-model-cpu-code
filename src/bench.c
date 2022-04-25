#include <assert.h>
#include <ctype.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "bench.h"
#include "bench_util.h"
#include "rng.h"
#include "schemes.h"

void check_function_pointer_exists(void *step_func, void *multistep_func, int use_multistep, char *scheme_str)
{
    if (use_multistep) {
        if (multistep_func == NULL) {
            printf("ERROR: Model has no multistep function for scheme %s\n", scheme_str);
            exit(EXIT_FAILURE);
        }
    } else {
        if (step_func == NULL) {
            printf("ERROR: Model has no step function for scheme %s\n", scheme_str);
            exit(EXIT_FAILURE);
        }
    }
}

void bench_singlerep(const struct cellmodel *model, struct options *opt)
{
    int num_threads = detect_num_threads();
    printf("Using %d threads\n", num_threads);

    unsigned long num_cells = opt->num_cells;
    int num_steps = opt->num_timesteps;
    double dt = opt->dt;
    int use_multistep = opt->use_multistep;

    unsigned int num_states = model->num_states;
    unsigned int num_parameters = model->num_parameters;

    size_t alignment_bytes = CELLMODEL_STATES_ALIGNMENT_BYTES;
    unsigned int padded_num_cells = (unsigned int) ceil_to_multiple_uint64(
            num_cells, alignment_bytes / sizeof(cellmodel_float_t));

    size_t states_size = padded_num_cells * num_states * sizeof(cellmodel_float_t);
    size_t parameters_2d_size = padded_num_cells * num_parameters * sizeof(cellmodel_float_t);
    size_t parameters_size = num_parameters * sizeof(cellmodel_float_t);
    size_t total_size = states_size + parameters_size;

    printf("Memory footprint: %lu bytes\n", total_size);

    cellmodel_float_t *states = aligned_alloc(alignment_bytes, states_size);
    //cellmodel_float_t *parameters = aligned_alloc(alignment_bytes, parameters_size);
    cellmodel_float_t *states_tmp1 = NULL;
    cellmodel_float_t *states_tmp2 = NULL;

    cellmodel_float_t *parameters = malloc(parameters_size);
    cellmodel_float_t *parameters_2d = NULL;

    if (use_multistep) {
        parameters_2d = aligned_alloc(alignment_bytes, parameters_2d_size);
    }

    model->init_states(states, num_cells, padded_num_cells);
    model->init_parameters(parameters);

    int is_stimulated_index = model->parameter_index("is_stimulated");
    if (is_stimulated_index == -1) {
        printf("Warning: Cannot apply stimulus when model has no is_stimulated parameter\n");
    }

    if (use_multistep) {
        assert(parameters_2d != NULL);
        for (int p = 0; p < num_parameters; p++) {
            #pragma omp parallel for
            for (int i = 0; i < num_cells; i++) {
                parameters_2d[p * padded_num_cells + i] = parameters[p];
            }
        }
    }

    const int num_V_names = 2;
    int V_index = -1;
    const char V_names[][10] = {"V_m", "V"};
    for (int i = 0; i < num_V_names; i++) {
        V_index = model->state_index(V_names[i]);
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
        for (int i = 0; i < num_cells; i++) {
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

    int stim_start_index = model->parameter_index("stim_start");
    if (stim_start_index != -1) {
        parameters[stim_start_index] = 1;
    }

    struct timespec timestamp_pre, timestamp_post;

    double t = 0;
    clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_pre);
    step_func step = NULL;
    multistep_func multistep = NULL;
    switch (opt->scheme) {
    case SCHEME_FE:
        step = model->step_FE;
        multistep = model->multistep_FE;
        check_function_pointer_exists(step, multistep, use_multistep, "FE");
        break;
    case SCHEME_RL:
        step = model->step_RL;
        multistep = model->multistep_RL;
        check_function_pointer_exists(step, multistep, use_multistep, "RL");
        break;
    case SCHEME_GRL1:
        step = model->step_GRL1;
        multistep = model->multistep_GRL1;
        check_function_pointer_exists(step, multistep, use_multistep, "GRL1");
        break;
    default:
        printf("ERROR: Unknown scheme with value %d\n", opt->scheme);
        exit(EXIT_FAILURE);
        break;
    }

    if (use_multistep) {
        assert(multistep != NULL);
        multistep(states, t, dt, parameters_2d, num_cells, padded_num_cells, num_steps);
    } else if (step) {
        // FE or GRL1
        for (int i = 0; i < num_steps; i++) {
            if (is_stimulated_index != -1) {
                if (t >= 1 && t < 2) {
                    parameters[is_stimulated_index] = 1;
                } else {
                    parameters[is_stimulated_index] = 0;
                }
            }
            step(states, t, dt, parameters, num_cells, padded_num_cells);
            t += dt;
        }
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

    double V_val = 0.0;
    if (model->layout == LAYOUT_STRUCT_OF_ARRAYS) {
        V_val = states[padded_num_cells * V_index];
    } else if (model->layout == LAYOUT_ARRAY_OF_STRUCTS) {
        V_val = states[V_index];
    }
    printf("V = %.15g at t = %g\n", V_val, t);

    free(states);
    free(parameters);
    if (use_multistep) {
        free(parameters_2d);
    }
    if (states_tmp1) {
        free(states_tmp1);
    }
    if (states_tmp2) {
        free(states_tmp2);
    }
}

void bench_reps(const struct cellmodel *model, struct options *opt)
{
    int num_threads = detect_num_threads();
    printf("Using %d threads\n", num_threads);

    unsigned long num_cells = opt->num_cells;
    int num_steps = opt->num_timesteps;
    double dt = opt->dt;
    int repetitions = opt->num_repetitions;
    int use_multistep = opt->use_multistep;

    unsigned int num_states = model->num_states;
    unsigned int num_parameters = model->num_parameters;

    size_t alignment_bytes = CELLMODEL_STATES_ALIGNMENT_BYTES;
    unsigned int padded_num_cells = (unsigned int) ceil_to_multiple_uint64(
            num_cells, alignment_bytes / sizeof(cellmodel_float_t));

    size_t states_size = padded_num_cells * num_states * sizeof(cellmodel_float_t);
    size_t parameters_2d_size = padded_num_cells * num_parameters * sizeof(cellmodel_float_t);
    size_t parameters_size = num_parameters * sizeof(cellmodel_float_t);
    size_t total_size = states_size + parameters_size;

    printf("Memory footprint: %lu bytes\n", total_size);


    int stim_start_index = model->parameter_index("stim_start");

    int is_stimulated_index = model->parameter_index("is_stimulated");
    if (stim_start_index == -1 && is_stimulated_index == -1) {
        printf("Warning: Cannot apply stimulus when model has no 'stim_start' or 'is_stimulated' parameter\n");
    }

    const int num_V_names = 2;
    int V_index = -1;
    const char V_names[][10] = {"V_m", "V"};
    for (int i = 0; i < num_V_names; i++) {
        V_index = model->state_index(V_names[i]);
        if (V_index != -1) {
            break;
        }
    }

    assert(V_index != -1);

    cellmodel_float_t *states = aligned_alloc(alignment_bytes, states_size);
    //cellmodel_float_t *parameters = aligned_alloc(alignment_bytes, parameters_size);
    cellmodel_float_t *states_tmp1 = NULL;
    cellmodel_float_t *states_tmp2 = NULL;

    cellmodel_float_t *parameters_2d = NULL;

    if (use_multistep) {
        parameters_2d = aligned_alloc(alignment_bytes, parameters_2d_size);
    }

    // select step function
    step_func step = NULL;
    multistep_func multistep = NULL;
    switch (opt->scheme) {
    case SCHEME_FE:
        step = model->step_FE;
        multistep = model->multistep_FE;
        check_function_pointer_exists(step, multistep, use_multistep, "FE");
        break;
    case SCHEME_RL:
        step = model->step_RL;
        multistep = model->multistep_RL;
        check_function_pointer_exists(step, multistep, use_multistep, "RL");
        break;
    case SCHEME_GRL1:
        step = model->step_GRL1;
        multistep = model->multistep_GRL1;
        check_function_pointer_exists(step, multistep, use_multistep, "GRL1");
        break;
    default:
        printf("ERROR: Unknown scheme with value %d\n", opt->scheme);
        exit(EXIT_FAILURE);
        break;
    }

    cellmodel_float_t *parameters = malloc(parameters_size);

    double *time_elapsed_arr = malloc(sizeof(double) * repetitions);
    double *cellsteps_per_second_arr = malloc(sizeof(double) * repetitions);
    double *mu_min_arr = malloc(sizeof(double) * repetitions);

    double t = -1;

    for (int r = 0; r < repetitions; r++) {
        model->init_states(states, num_cells, padded_num_cells);
        model->init_parameters(parameters);

        if (use_multistep) {
            assert(parameters_2d != NULL);
            for (int p = 0; p < num_parameters; p++) {
                #pragma omp parallel for
                for (int i = 0; i < num_cells; i++) {
                    parameters_2d[p * padded_num_cells + i] = parameters[p];
                }
            }
        }

        if (stim_start_index != -1) {
            parameters[stim_start_index] = 1;
        }

        // initialise V
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
            for (int i = 0; i < num_cells; i++) {
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

        t = 0;
        struct timespec timestamp_pre, timestamp_post;
        clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_pre);
        if (use_multistep) {
            assert(multistep != NULL);
            multistep(states, t, dt, parameters_2d, num_cells, padded_num_cells, num_steps);
        } else if (step) {
            // FE, RL or GRL1
            for (int i = 0; i < num_steps; i++) {
                if (is_stimulated_index != -1) {
                    if (t >= 1 && t < 2) {
                        parameters[is_stimulated_index] = 1;
                    } else {
                        parameters[is_stimulated_index] = 0;
                    }
                }
                step(states, t, dt, parameters, num_cells, padded_num_cells);
                t += dt;
            }
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

    double V_val = 0.0;
    if (model->layout == LAYOUT_STRUCT_OF_ARRAYS) {
        V_val = states[padded_num_cells * V_index];
    } else if (model->layout == LAYOUT_ARRAY_OF_STRUCTS) {
        V_val = states[V_index];
    }
    printf("V = %.15g at t = %g\n", V_val, t);

    free(states);
    free(parameters);
    if (use_multistep) {
        free(parameters_2d);
    }
    if (states_tmp1) {
        free(states_tmp1);
    }
    if (states_tmp2) {
        free(states_tmp2);
    }

    free(time_elapsed_arr);
    free(cellsteps_per_second_arr);
    free(mu_min_arr);
}

void bench(const struct cellmodel *model, struct options *opt)
{
    if (opt->num_repetitions > 1) {
        bench_reps(model, opt);
    } else {
        bench_singlerep(model, opt);
    }
}
