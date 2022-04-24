
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "TP06_naive.h"
#include "TP06_simd.h"
#include "common.h"
#include "schemes.h"
#include "solve.h"
#include "trace.h"
#include "../../dump_trace_cmdline.h"


static int parse_arguments_ggo(int argc, char **argv, struct gengetopt_args_info *args_info_ptr)
{
    int result = 0;

    struct cmdline_parser_params *params;


    /* initialize the parameters structure */
    params = cmdline_parser_params_create();

    /* call the command line parser */
    if (cmdline_parser(argc, argv, args_info_ptr) != 0) {
        result = 1;
        goto stop;
    }

    // the argument structure has already been initialised when parsing the command line arguments
    params->initialize = 0;
    // any options provided on the command line should take priority over the config file
    params->override = 0;

    free(params);
    return result;

stop:
    /* deallocate structures */
    cmdline_parser_free(args_info_ptr);
    free(params);
    return result;
}

int main(int argc, char *argv[])
{
    struct gengetopt_args_info *args_info = malloc(sizeof(struct gengetopt_args_info));
    if (parse_arguments_ggo(argc, argv, args_info)) {
        printf("ERROR: Failed to parse command line arguments\n");
        exit(EXIT_FAILURE);
    }

    const struct cellmodel *model = NULL;
    double T_end = 1000;
    double solve_dt = 1E-3;
    int store_period = 1;
    scheme_type scheme = SCHEME_RL;
    int imp_str_size = 20;
    char imp_str[imp_str_size];
    char *output_directory = args_info->output_directory_arg;
    create_dir_if_missing(output_directory);

    switch (args_info->implementation_arg) {
        case implementation__NULL:
        case implementation_arg_naive:
            model = &model_TP06_naive;
            snprintf(imp_str, imp_str_size, "naive");
            break;
        case implementation_arg_simd:
            model = &model_TP06_simd;
            snprintf(imp_str, imp_str_size, "simd");
            break;
    }
    switch (args_info->scheme_arg) {
        case scheme__NULL:
        case scheme_arg_FE:
            scheme = SCHEME_FE;
            break;
        case scheme_arg_RL:
            scheme = SCHEME_RL;
            break;
        case scheme_arg_GRL1:
            scheme = SCHEME_GRL1;
            break;

    }

    if (args_info->solve_dt_given) {
        solve_dt = args_info->solve_dt_arg;
    }
    if (args_info->store_period_given) {
        store_period = args_info->store_period_arg;
    }


    cellmodel_float_t *V_trace;
    cellmodel_float_t *t_trace;
    int trace_length;

    long num_timesteps = lround(T_end / solve_dt);

    printf("Solving with %ld steps\n", num_timesteps);

    solve_single(model, solve_dt, T_end, store_period, scheme, &V_trace, &t_trace, &trace_length);

    printf("Storing trace of length %d\n", trace_length);

    size_t scheme_str_size = 20;
    char scheme_str[scheme_str_size];
    get_scheme_str(scheme_str, scheme_str_size, scheme);

    size_t basename_size = 1024;
    char subdir[basename_size];
    char basename[basename_size];
    snprintf(subdir, basename_size, "%s/%s", output_directory, imp_str);
    create_dir_if_missing(subdir);
    snprintf(basename, basename_size, "%s/%s/TP06_%s_%g", output_directory, imp_str, scheme_str, solve_dt);

    printf("Storing trace with basename '%s'\n", basename);

    dump_trace(t_trace, V_trace, trace_length, solve_dt, scheme, basename);

    cmdline_parser_free(args_info);
    free(args_info);
}
