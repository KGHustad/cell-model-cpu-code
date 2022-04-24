#include "bench_util.h"

#include <ctype.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rng.h"

static const char *scheme_opts[] = {"fe", "rl", "grl1", "grl2", NULL};
static const char *v_var_opts[] = {"homogeneous", "noise", "linear", NULL};

static void print_usage(int argc, char *argv[])
{
    printf("Usage: %s \n"
           "    [-s {FE,GRL1}] [--scheme {FE,GRL1}]\n"
           "    [-t <int>] [--timesteps <int>]\n"
           "    [-N <int>] [--num_cells <int>]\n"
           "    [-d <float>] [--dt <float>]\n"
           "    [-v {homogeneous,noise,linear}] [--v_variation {...}]\n"
           "    [--seed <unsigned int>]\n"
           "    [-r <int>] [--repetitions <int>]\n",
           argv[0]);
}

enum {
    V_STEP_SYM = 256,
    CA_STEP_SYM = 257,
    SEED_SYM = 258,
};

void bench_parse_args(int argc, char *argv[], struct options *opt)
{
    opt->verbose_flag = 1;
    opt->num_timesteps = 5;
    opt->num_repetitions = 1;
    opt->scheme = SCHEME_FE;
    opt->num_cells = (unsigned int) 11.5E6;
    opt->dt = 1E-4;
    opt->use_multistep = 0;
    opt->v_var = V_HOMOGENEOUS;
    opt->seed = rng_draw_seed();
    int c;
    while (1) {

        struct option long_options[] = {/* These options set a flag. */
                                        {"verbose", no_argument, &opt->verbose_flag, 1},
                                        {"quiet", no_argument, &opt->verbose_flag, 0},
                                        {"multistep", no_argument, &opt->use_multistep, 1},
                                        /* These options don’t set a flag.
                We distinguish them by their indices. */
                                        {"help", required_argument, 0, 'h'},
                                        {"dt", required_argument, 0, 'd'},
                                        {"scheme", required_argument, 0, 's'},
                                        {"timesteps", required_argument, 0, 't'},
                                        {"num_cells", required_argument, 0, 'N'},
                                        {"v_variation", required_argument, 0, 'v'},
                                        {"seed", required_argument, 0, SEED_SYM},
                                        {"V_step", required_argument, 0, V_STEP_SYM},
                                        {"Ca_step", required_argument, 0, CA_STEP_SYM},
                                        {"repetitions", required_argument, 0, 'r'},
                                        {0, 0, 0, 0}};

        int option_index = 0;
        c = getopt_long(argc, argv, "hd:r:s:t:N:", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c) {
        case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
                break;
            printf("option %s", long_options[option_index].name);
            if (optarg)
                printf(" with arg %s", optarg);
            printf("\n");
            break;

        case 'h':
            print_usage(argc, argv);
            exit(EXIT_SUCCESS);
            break;

        case 'd': {
            char *endp;
            double val = strtod(optarg, &endp);
            if (endp == optarg) {
                printf("error: dt could not be interpreted as number, got '%s'\n", optarg);
                exit(EXIT_FAILURE);
            }
            if (val <= 0) {
                printf("error: dt must be positive, got '%s'\n", optarg);
                exit(EXIT_FAILURE);
            }
            opt->dt = val;
            break;
        }

        case 's': {
            int k;
            size_t optarg_size = strlen(optarg) + 1;
            char optarg_lowercase[optarg_size];
            strcpy(optarg_lowercase, optarg);
            for (int i = 0; i < optarg_size - 1; i++) {
                optarg_lowercase[i] = tolower(optarg[i]);
            }
            for (k = 0; k < NUM_SCHEMES; k++) {
                if (strcmp(scheme_opts[k], optarg_lowercase) == 0) {
                    break;
                }
            }
            if (k == NUM_SCHEMES) {
                printf("error: Unrecognised value for --scheme option: %s\n", optarg);
                exit(EXIT_FAILURE);
            }
            opt->scheme = (scheme_type) k;
            break;
        }

        case 'r': {
            int val = atoi(optarg);
            if (val <= 0) {
                printf("error: number of repetitions must be positive, got '%s'\n", optarg);
                exit(EXIT_FAILURE);
            }
            opt->num_repetitions = val;
            break;
        }

        case 't': {
            int val = atoi(optarg);
            if (val <= 0) {
                printf("error: number of timesteps must be positive, got '%s'\n", optarg);
                exit(EXIT_FAILURE);
            }
            opt->num_timesteps = val;
            break;
        }

        case 'N': {
            char *endp = optarg;
            double double_val = strtod(optarg, &endp);
            if (endp == optarg) {
                printf("error: num_cells could not be interpreted as number, got '%s'\n", optarg);
                exit(EXIT_FAILURE);
            }
            unsigned int val = (unsigned int) (double_val + 0.5);
            if (val <= 0) {
                printf("error: num_cells must be positive, got '%s'\n", optarg);
                exit(EXIT_FAILURE);
            }
            opt->num_cells = val;
            break;
        }

        case 'v': {
            int k;
            size_t optarg_size = strlen(optarg) + 1;
            char optarg_lowercase[optarg_size];
            strcpy(optarg_lowercase, optarg);
            for (int i = 0; i < optarg_size - 1; i++) {
                optarg_lowercase[i] = tolower(optarg[i]);
            }
            for (k = 0; k < NUM_V_VAR; k++) {
                if (strcmp(v_var_opts[k], optarg_lowercase) == 0) {
                    break;
                }
            }
            if (k == NUM_V_VAR) {
                printf("error: Unrecognised value for --v_variation option: %s\n", optarg);
                exit(EXIT_FAILURE);
            }
            opt->v_var = (enum v_variation) k;
            break;
        }

        case SEED_SYM: {
            char *endp = optarg;
            unsigned long val = strtoul(optarg, &endp, 10);
            if (endp == optarg) {
                printf("error: seed could not be interpreted as unsigned integer, got '%s'\n",
                       optarg);
                exit(EXIT_FAILURE);
            }
            opt->seed = val;
            break;
        }

        case V_STEP_SYM: {
            char *endp;
            double val = strtod(optarg, &endp);
            if (endp == optarg) {
                printf("error: V_step could not be interpreted as number, got '%s'\n", optarg);
                exit(EXIT_FAILURE);
            }
            if (val <= 0) {
                printf("error: V_step must be positive, got '%s'\n", optarg);
                exit(EXIT_FAILURE);
            }
            opt->V_step = val;
            break;
        }

        case CA_STEP_SYM: {
            char *endp;
            double val = strtod(optarg, &endp);
            if (endp == optarg) {
                printf("error: Ca_step could not be interpreted as number, got '%s'\n", optarg);
                exit(EXIT_FAILURE);
            }
            if (val <= 0) {
                printf("error: Ca_step must be positive, got '%s'\n", optarg);
                exit(EXIT_FAILURE);
            }
            opt->Ca_step = val;
            break;
        }

        case '?':
            /* getopt_long already printed an error message. */
            exit(EXIT_FAILURE);
            break;
        }
    }

    if (opt->verbose_flag) {
        printf("Number of cells: %lu\n", opt->num_cells);
        printf("dt: %g\n", opt->dt);
        switch (opt->scheme) {
        case SCHEME_FE:
            printf("Scheme: FE\n");
            break;
        case SCHEME_GRL1:
            printf("Scheme: GRL1\n");
            break;
        default:
            break;
        }
    }
}

uint64_t ceil_div_uint64(uint64_t a, uint64_t b)
{
    return (a + b - 1) / b;
}

uint64_t ceil_to_multiple_uint64(uint64_t a, uint64_t b)
{
    return ceil_div_uint64(a, b) * b;
}

int detect_num_threads()
{
    int count = 0;
    #pragma omp parallel
    {
        #pragma omp critical
        {
            count++;
        }
    }
    return count;
}

int compare_double(const void *a_void, const void *b_void)
{
    const double a = *(const double *) a_void;
    const double b = *(const double *) b_void;
    if (a > b) {
        return 1;
    } else if (a < b) {
        return -1;
    } else {
        return 0;
    }
}

int compare_double_inv(const void *a_void, const void *b_void)
{
    return compare_double(b_void, a_void);
}

double relative_diff(double a, double b)
{
    return (a - b) / b;
}
