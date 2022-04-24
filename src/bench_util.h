#ifndef _BENCH_UTIL_H
#define _BENCH_UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

#include "cellmodel.h"

enum v_variation {
    V_HOMOGENEOUS,
    V_NOISE,
    V_LINEAR,
    NUM_V_VAR,
};

struct options {
    int verbose_flag;
    int num_timesteps;
    int num_repetitions;
    scheme_type scheme;
    unsigned long num_cells;
    int use_multistep;
    double dt;
    enum v_variation v_var; // variation between cells in the initial values of V
    uint64_t seed;
    double V_step;
    double Ca_step;
};

void bench_parse_args(int argc, char *argv[], struct options *opt);

uint64_t ceil_div_uint64(uint64_t a, uint64_t b);
uint64_t ceil_to_multiple_uint64(uint64_t a, uint64_t b);
int detect_num_threads();
int compare_double(const void *a_void, const void *b_void);
int compare_double_inv(const void *a_void, const void *b_void);
double relative_diff(double a, double b);

#ifdef __cplusplus
} // extern "C"
#endif

#endif
