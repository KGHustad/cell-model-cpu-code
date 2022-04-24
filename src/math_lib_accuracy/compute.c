#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "common.h"

void apply_exp(const double *__restrict in, double *__restrict out, int N)
{
#ifdef VECTOR_LENGTH
    #pragma omp for simd simdlen(VECTOR_LENGTH)
#else
    #pragma omp for simd
#endif
    for (int i = 0; i < N; i++) {
        out[i] = exp(in[i]);
    }
}

void apply_expm1(const double *__restrict in, double *__restrict out, int N)
{
#ifdef VECTOR_LENGTH
    #pragma omp for simd simdlen(VECTOR_LENGTH)
#else
    #pragma omp for simd
#endif
    for (int i = 0; i < N; i++) {
        out[i] = expm1(in[i]);
    }
}

void apply_log(const double *__restrict in, double *__restrict out, int N)
{
#ifdef VECTOR_LENGTH
    #pragma omp for simd simdlen(VECTOR_LENGTH)
#else
    #pragma omp for simd
#endif
    for (int i = 0; i < N; i++) {
        out[i] = log(in[i]);
    }
}

void apply_pow(const double *__restrict in_base, const double *__restrict in_exponent,
               double *__restrict out, int N)
{
#ifdef VECTOR_LENGTH
    #pragma omp for simd simdlen(VECTOR_LENGTH)
#else
    #pragma omp for simd
#endif
    for (int i = 0; i < N; i++) {
        out[i] = pow(in_base[i], in_exponent[i]);
    }
}

int dump_array_to_binary_file(char *filename, double *b, int N)
{
    assert(b != NULL);

    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        fprintf(stderr, "Error: Could not open binary file at '%s' for writing\n", filename);
        exit(EXIT_FAILURE);
    }

    size_t elements_written = fwrite(b, sizeof(double), N, f);
    assert(elements_written == N);

    fclose(f);

    return 0;
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <basename>\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    char *basename = argv[1];

    // load init file with array length
    int N;
    load_init_file(basename, &N);
    printf("N=%d\n", N);

    // allocate arrays
    size_t array_size = sizeof(double) * N;

    double *exp_input = malloc(array_size);
    double *exp_output = malloc(array_size);

    double *expm1_input = malloc(array_size);
    double *expm1_output = malloc(array_size);

    double *log_input = malloc(array_size);
    double *log_output = malloc(array_size);

    double *pow_input_b = malloc(array_size);
    double *pow_input_e = malloc(array_size);
    double *pow_output = malloc(array_size);

    // load array contents from disk
    int filename_size = 1024;
    char filename[filename_size];

    snprintf(filename, filename_size, "%s_%s.bin", basename, "in_exp");
    load_binary_file_into_array(filename, exp_input, N);

    snprintf(filename, filename_size, "%s_%s.bin", basename, "in_expm1");
    load_binary_file_into_array(filename, expm1_input, N);

    snprintf(filename, filename_size, "%s_%s.bin", basename, "in_log");
    load_binary_file_into_array(filename, log_input, N);

    snprintf(filename, filename_size, "%s_%s.bin", basename, "in_pow_base");
    load_binary_file_into_array(filename, pow_input_b, N);

    snprintf(filename, filename_size, "%s_%s.bin", basename, "in_pow_exponent");
    load_binary_file_into_array(filename, pow_input_e, N);

    // compute
    apply_exp(exp_input, exp_output, N);
    apply_expm1(expm1_input, expm1_output, N);
    apply_log(log_input, log_output, N);
    apply_pow(pow_input_b, pow_input_e, pow_output, N);

    // dump output arrays to disk
    snprintf(filename, filename_size, "%s_%s.bin", basename, "out_exp");
    dump_array_to_binary_file(filename, exp_output, N);

    snprintf(filename, filename_size, "%s_%s.bin", basename, "out_expm1");
    dump_array_to_binary_file(filename, expm1_output, N);

    snprintf(filename, filename_size, "%s_%s.bin", basename, "out_log");
    dump_array_to_binary_file(filename, log_output, N);

    snprintf(filename, filename_size, "%s_%s.bin", basename, "out_pow");
    dump_array_to_binary_file(filename, pow_output, N);


    // free arrays
    free(exp_input);
    free(exp_output);

    free(expm1_input);
    free(expm1_output);

    free(log_input);
    free(log_output);

    free(pow_input_b);
    free(pow_input_e);
    free(pow_output);

    return EXIT_SUCCESS;
}
