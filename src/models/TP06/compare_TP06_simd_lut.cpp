#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../compare_naive_to_lut.hpp"

#include "TP06_naive.h"
#include "TP06_simd_lut.hpp"

struct lut_options {
    double V_step;
    double Ca_step;
};

static void print_usage(int argc, char *argv[])
{
    printf("Usage: %s \n"
           "    [-V <float] [--V_step <float>]\n"
           "    [-C <float] [--Ca_step <float>]\n",
           argv[0]);
}

void parse_options(int argc, char *argv[], struct lut_options *opt)
{
    int c;
    while (1) {
        struct option long_options[] = {/* These options set a flag. */
                                        /* These options don’t set a flag.
                We distinguish them by their indices. */
                                        {"help", required_argument, 0, 'h'},
                                        {"V_step", required_argument, 0, 'V'},
                                        {"Ca_step", required_argument, 0, 'C'},
                                        {0, 0, 0, 0}};

        int option_index = 0;
        c = getopt_long(argc, argv, "hV:C:", long_options, &option_index);

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

        case 'V': {
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

        case 'C': {
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
}

int main(int argc, char *argv[])
{
    long num_cells = (long) 1;
    int num_steps = 800000;

    struct lut_options opt;

    opt.V_step = 0.05;
    opt.Ca_step = 0.00001;

    parse_options(argc, argv, &opt);

    double V_step = opt.V_step;
    double Ca_step = opt.Ca_step;

    printf("LUT method: Linear interpolation\n");

    printf("V_step: %g\n", V_step);
    printf("Ca_step: %g\n", Ca_step);

    compare_naive_to_lut(&model_TP06_naive, &model_TP06_simd_lut, num_cells, num_steps, V_step, Ca_step,
                         SCHEME_RL);
}
