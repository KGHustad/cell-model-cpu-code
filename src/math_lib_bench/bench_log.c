#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "common.h"

void init_log(int N, double *__restrict a, double *__restrict b)
{
    for (int i = 0; i < N; i++) {
        a[i] = 1 + ((double) i) / N;
        b[i] = 0;
    }
}

#pragma GCC optimize("-fno-tree-vectorize")
void perturb_scalar(int N, double *__restrict b)
{
    double eps = 1E-6;
    #pragma clang loop vectorize(disable)
    #pragma novector
    for (int i = 0; i < N; i++) {
        b[i] += eps;
    }
}
#pragma GCC reset_options

void perturb_omp_simd(int N, double *__restrict b)
{
    double eps = 1E-6;
#ifdef VECTOR_LENGTH
    #pragma omp for simd simdlen(VECTOR_LENGTH)
#else
    #pragma omp for simd
#endif
    for (int i = 0; i < N; i++) {
        b[i] += eps;
    }
}

#pragma GCC optimize("-fno-tree-vectorize")
void invert_scalar(int N, double *__restrict a, double *__restrict b)
{
    /*
    y = log(x)
    x = exp(y)
    */
    #pragma clang loop vectorize(disable)
    #pragma novector
    for (int i = 0; i < N; i++) {
        a[i] = exp(b[i]);
    }
}
#pragma GCC reset_options

void invert_omp_simd(int N, double *__restrict a, double *__restrict b)
{
    /*
    y = log(x)
    x = exp(y)
    */
#ifdef VECTOR_LENGTH
    #pragma omp for simd simdlen(VECTOR_LENGTH)
#else
    #pragma omp for simd
#endif
    for (int i = 0; i < N; i++) {
        a[i] = exp(b[i]);
    }
}

#pragma GCC optimize("-fno-tree-vectorize")
void bench_log_scalar(int N, double *__restrict a, double *__restrict b)
{
    #pragma clang loop vectorize(disable)
    #pragma novector
    for (int i = 0; i < N; i++) {
        b[i] = log(a[i]);
    }
}
#pragma GCC reset_options

double array_mean(int N, double *__restrict a)
{
    double mean = 0;
    double N_inv = 1. / N;
    for (int i = 0; i < N; i++) {
        mean += N_inv * a[i];
    }
    return mean;
}

void bench_log_omp_simd(int N, double *__restrict a, double *__restrict b)
{
#ifdef VECTOR_LENGTH
    #pragma omp for simd simdlen(VECTOR_LENGTH)
#else
    #pragma omp for simd
#endif
    for (int i = 0; i < N; i++) {
        b[i] = log(a[i]);
    }
}

#ifdef HAS_ACCELERATE_FRAMEWORK
void bench_log_vforce(int N, double *__restrict a, double *__restrict b)
{
    vvlog(b, a, &N);
}
#endif

void bench_log(int N)
{
    printf("Benchmarking math function log(). Array length: %d\n", N);

    double *a = malloc(sizeof(double) * N);
    double *b = malloc(sizeof(double) * N);

    double time_elapsed, throughput, mu;
    double mean;

    struct timespec timestamp_pre, timestamp_post;

    // run scalar code
    init_log(N, a, b);
    clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_pre);
    bench_log_scalar(N, a, b);
    clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_post);
    mean = array_mean(N, b);
    time_elapsed = timestamp_post.tv_sec - timestamp_pre.tv_sec + 1E-9*(timestamp_post.tv_nsec - timestamp_pre.tv_nsec);
    throughput = N / time_elapsed;
    mu = 2*sizeof(double) * throughput;

    printf("Variant: Scalar. Mean: %g  Time elapsed: %g seconds. Throughput: %g. Memory traffic: %.4g GB/s\n", mean, time_elapsed, throughput, 1E-9*mu);
    double throughput_scalar = throughput;

    // run simd code
    init_log(N, a, b);
    clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_pre);
    bench_log_omp_simd(N, a, b);
    clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_post);
    mean = array_mean(N, b);
    time_elapsed = timestamp_post.tv_sec - timestamp_pre.tv_sec + 1E-9*(timestamp_post.tv_nsec - timestamp_pre.tv_nsec);
    throughput = N / time_elapsed;
    mu = 2*sizeof(double) * throughput;
    printf("Variant: SIMD.   Mean: %g  Time elapsed: %g seconds. Throughput: %g. Memory traffic: %.4g GB/s\n", mean, time_elapsed, throughput, 1E-9*mu);
    double throughput_simd = throughput;

#ifdef HAS_ACCELERATE_FRAMEWORK
    // run vforce code
    init_log(N, a, b);
    clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_pre);
    bench_log_vforce(N, a, b);
    clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_post);
    mean = array_mean(N, b);
    time_elapsed = timestamp_post.tv_sec - timestamp_pre.tv_sec + 1E-9*(timestamp_post.tv_nsec - timestamp_pre.tv_nsec);
    throughput = N / time_elapsed;
    mu = 2*sizeof(double) * throughput;
    printf("Variant: vForce. Mean: %g  Time elapsed: %g seconds. Throughput: %g. Memory traffic: %.4g GB/s\n", mean, time_elapsed, throughput, 1E-9*mu);
    double throughput_vforce = throughput;
#endif

    printf("Speedup from scalar to SIMD: %g\n", throughput_simd / throughput_scalar);
#ifdef HAS_ACCELERATE_FRAMEWORK
    printf("Speedup from scalar to vForce: %g\n", throughput_vforce / throughput_scalar);
#endif

    free(a);
    free(b);
}

void bench_log_rep(int N, int repetitions)
{
    printf("Benchmarking math function log(). Array length: %d. Number of repetitions: %d\n", N, repetitions);

    double *a = malloc(sizeof(double) * N);
    double *b = malloc(sizeof(double) * N);

    double time_elapsed, throughput, mu;
    double mean;

    struct timespec timestamp_pre, timestamp_post;
    long second_counter;
    long nanosecond_counter;

    // run scalar code
    init_log(N, a, b);
    second_counter = 0;
    nanosecond_counter = 0;
    for (int r = 0; r < repetitions; r++) {
        clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_pre);
        bench_log_scalar(N, a, b);
        clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_post);
        second_counter += timestamp_post.tv_sec - timestamp_pre.tv_sec;
        nanosecond_counter += timestamp_post.tv_nsec - timestamp_pre.tv_nsec;

        // prepare for next rep
        perturb_scalar(N, b);
        invert_scalar(N, a, b);
    }
    mean = array_mean(N, b);
    time_elapsed = second_counter + 1E-9*(nanosecond_counter);
    throughput = N * (long) repetitions / time_elapsed;
    mu = 2*sizeof(double) * throughput;

    printf("Variant: Scalar. Mean: %g  Time elapsed: %g seconds. Throughput: %g. Memory traffic: %.4g GB/s\n", mean, time_elapsed, throughput, 1E-9*mu);
    double throughput_scalar = throughput;

    // run simd code
    init_log(N, a, b);
    second_counter = 0;
    nanosecond_counter = 0;
    for (int r = 0; r < repetitions; r++) {
        clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_pre);
        bench_log_omp_simd(N, a, b);
        clock_gettime(CLOCK_MONOTONIC_RAW, &timestamp_post);
        second_counter += timestamp_post.tv_sec - timestamp_pre.tv_sec;
        nanosecond_counter += timestamp_post.tv_nsec - timestamp_pre.tv_nsec;

        // prepare for next rep
        perturb_omp_simd(N, b);
        invert_omp_simd(N, a, b);
    }
    mean = array_mean(N, b);
    time_elapsed = second_counter + 1E-9*(nanosecond_counter);
    throughput = N * (long) repetitions / time_elapsed;
    mu = 2*sizeof(double) * throughput;
    printf("Variant: SIMD.   Mean: %g  Time elapsed: %g seconds. Throughput: %g. Memory traffic: %.4g GB/s\n", mean, time_elapsed, throughput, 1E-9*mu);
    double throughput_simd = throughput;

    printf("Speedup from scalar to SIMD: %g\n", throughput_simd / throughput_scalar);

    free(a);
    free(b);
}

int main(int argc, char *argv[])
{
    int N = (int) 1E8;
    int repetitions = 1;
    if (argc > 1) {
        N = atoi(argv[1]);
    }
    if (argc > 2) {
        repetitions = atoi(argv[2]);
        assert(repetitions > 0);
    }

    if (repetitions > 1) {
        bench_log_rep(N, repetitions);
    } else {
        bench_log(N);
    }

    return EXIT_SUCCESS;
}
