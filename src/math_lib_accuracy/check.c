#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpfr.h>

#include "common.h"

typedef int (*univariate_mpfr_func)(mpfr_t, const mpfr_t, mpfr_rnd_t);
typedef int (*bivariate_mpfr_func)(mpfr_t, const mpfr_t, const mpfr_t, mpfr_rnd_t);

static const int mpfr_bits = 120;
static const mpfr_rnd_t rounding_mode = MPFR_RNDN;

int compare_double(const void *a_ptr, const void *b_ptr)
{
    const double a = *((const double *) a_ptr);
    const double b = *((const double *) b_ptr);
    if (a < b) {
        return -1;
    } else if (a > b) {
        return 1;
    } else {
        return 0;
    }
}


void check_univariate_f(const double *__restrict f_input, const double *__restrict f_output, int N,
                        univariate_mpfr_func f, char *function_name)
{
    printf("\nChecking %s()\n", function_name);

    mpfr_t m_input, m_output, m_rounding_error, m_approx_error, m_approx_error_ulps, m_ulps,
            m_nearest;
    mpfr_init2(m_input, mpfr_bits);
    mpfr_init2(m_output, mpfr_bits);
    mpfr_init2(m_rounding_error, mpfr_bits);
    mpfr_init2(m_approx_error, mpfr_bits);
    mpfr_init2(m_approx_error_ulps, mpfr_bits);
    mpfr_init2(m_ulps, mpfr_bits);
    mpfr_init2(m_nearest, mpfr_bits);

    double *approx_error_ulps_a = malloc(sizeof(double) * N);
    int equal_count = 0;
    double max_error_ulps = 0;
    for (int i = 0; i < N; i++) {
        double input_fp64 = f_input[i];
        double output_fp64_approx = f_output[i];

        // set input
        mpfr_set_d(m_input, input_fp64, rounding_mode);

        // evaluate function with MPFR
        //mpfr_exp(m_output, m_input, rounding_mode);
        (*f)(m_output, m_input, rounding_mode);

        // compute rounding error to find direction for ULPS
        double output_fp64_exact = mpfr_get_d(m_output, rounding_mode);
        if (output_fp64_approx == output_fp64_exact) {
            equal_count++;
        }
        mpfr_d_sub(m_rounding_error, output_fp64_exact, m_output, rounding_mode);

        // determine ULPS for output
        double rounding_error = mpfr_get_d(m_rounding_error, rounding_mode);
        double direction = rounding_error > 0 ? DBL_MAX : DBL_MIN;
        if (fabs(rounding_error) == 0) {
            direction = DBL_MAX;
        }
        double nearest_fp64 = nextafter(output_fp64_exact, direction);

        // compute approximation error
        mpfr_d_sub(m_approx_error, output_fp64_approx, m_output, rounding_mode);
        double approx_error_ulps;

        // ulps is subnormal, compute it with
        mpfr_set_d(m_nearest, nearest_fp64, rounding_mode);
        mpfr_sub_d(m_ulps, m_nearest, output_fp64_exact, rounding_mode);

        mpfr_div(m_approx_error_ulps, m_approx_error, m_ulps, rounding_mode);
        approx_error_ulps = mpfr_get_d(m_approx_error_ulps, rounding_mode);
        if (fabs(approx_error_ulps) > 4) {
            mpfr_dump(m_ulps);
            printf("i=%d  rounding_error=%g  input=%.17e  approx_error_ulps=%g\n"
                   "nearest_fp64 =       %.17e\n"
                   "output_fp64_exact =  %.17e\n"
                   "output_fp64_approx = %.17e\n\n",
                   i, rounding_error, input_fp64, approx_error_ulps, nearest_fp64,
                   output_fp64_exact, output_fp64_approx);
        }

        approx_error_ulps_a[i] = fabs(approx_error_ulps);
        if (fabs(approx_error_ulps) > max_error_ulps) {
            max_error_ulps = fabs(approx_error_ulps);
        }
    }
    int unequal_count = N - equal_count;
    double unequal_ratio = unequal_count / (double) N;
    printf("%7d / %7d (%g%%) were bitwise unequal\n", unequal_count, N, 100 * unequal_ratio);
    printf("Max approximation error: %g ULPS\n", max_error_ulps);

#if 0
    qsort(approx_error_ulps_a, N, sizeof(double), &compare_double);

    double error_percentile_25 = approx_error_ulps_a[(N*25L)/100];
    double error_percentile_50 = approx_error_ulps_a[(N*50L)/100];
    double error_percentile_75 = approx_error_ulps_a[(N*75L)/100];
    double error_percentile_99 = approx_error_ulps_a[(N*99L)/100];
    printf("25th percentile: %g ULPS\n", error_percentile_25);
    printf("50th percentile: %g ULPS\n", error_percentile_50);
    printf("75th percentile: %g ULPS\n", error_percentile_75);
    printf("99th percentile: %g ULPS\n", error_percentile_99);
    printf("Max: %g ULPS\n", approx_error_ulps_a[N-1]);

    for (int i = N-10; i < N; i++) {
        printf("approx_error_ulps_a[%d]: %g ULPS\n", i, approx_error_ulps_a[i]);
    }
#endif


    free(approx_error_ulps_a);

    mpfr_clear(m_input);
    mpfr_clear(m_output);
    mpfr_clear(m_rounding_error);
    mpfr_clear(m_approx_error);
    mpfr_clear(m_approx_error_ulps);
    mpfr_clear(m_ulps);
    mpfr_clear(m_nearest);
    mpfr_free_cache();
}

void check_pow(const double *__restrict f_input_base, const double *__restrict f_input_exponent,
               const double *__restrict f_output, int N)
{
    printf("\nChecking pow()\n");

    mpfr_t m_input_base, m_input_exponent, m_output, m_rounding_error, m_approx_error,
            m_approx_error_ulps, m_ulps, m_nearest;
    mpfr_init2(m_input_base, mpfr_bits);
    mpfr_init2(m_input_exponent, mpfr_bits);
    mpfr_init2(m_output, mpfr_bits);
    mpfr_init2(m_rounding_error, mpfr_bits);
    mpfr_init2(m_approx_error, mpfr_bits);
    mpfr_init2(m_approx_error_ulps, mpfr_bits);
    mpfr_init2(m_ulps, mpfr_bits);
    mpfr_init2(m_nearest, mpfr_bits);

    double *approx_error_ulps_a = malloc(sizeof(double) * N);
    int equal_count = 0;
    int both_zero_count = 0;
    int subnormal_answer_count = 0;
    double max_error_ulps = 0;
    for (int i = 0; i < N; i++) {
        double input_base_fp64 = f_input_base[i];
        double input_exponent_fp64 = f_input_exponent[i];
        double output_fp64_approx = f_output[i];

        // set input
        mpfr_set_d(m_input_base, input_base_fp64, rounding_mode);
        mpfr_set_d(m_input_exponent, input_exponent_fp64, rounding_mode);

        // evaluate function with MPFR
        mpfr_pow(m_output, m_input_base, m_input_exponent, rounding_mode);

        // compute rounding error to find direction for ULPS
        double output_fp64_exact = mpfr_get_d(m_output, rounding_mode);
        if (output_fp64_approx == output_fp64_exact) {
            equal_count++;
            if (output_fp64_exact == 0) {
                // it makes no sense to compute ULPS error when both values are zero
                approx_error_ulps_a[i] = 0;
                both_zero_count++;
                continue;
            }
        }

        if (fabs(output_fp64_exact) < 2.5E-308 && fabs(output_fp64_approx) == 0) {
            approx_error_ulps_a[i] = 0;
            subnormal_answer_count++;
            continue;
        }

        mpfr_d_sub(m_rounding_error, output_fp64_exact, m_output, rounding_mode);

        // determine ULPS for output
        double rounding_error = mpfr_get_d(m_rounding_error, rounding_mode);
        double direction = rounding_error > 0 ? DBL_MAX : DBL_MIN;
        if (fabs(rounding_error) == 0) {
            direction = DBL_MAX;
        }
        double nearest_fp64 = nextafter(output_fp64_exact, direction);

        // compute approximation error
        mpfr_d_sub(m_approx_error, output_fp64_approx, m_output, rounding_mode);

        mpfr_set_d(m_nearest, nearest_fp64, rounding_mode);
        mpfr_sub_d(m_ulps, m_nearest, output_fp64_exact, rounding_mode);
        if (mpfr_zero_p(m_ulps)) {
            approx_error_ulps_a[i] = 0;
            continue;
        }

        mpfr_div(m_approx_error_ulps, m_approx_error, m_ulps, rounding_mode);
        double approx_error_ulps = mpfr_get_d(m_approx_error_ulps, rounding_mode);
        if (fabs(approx_error_ulps) > 4) {
            //mpfr_dump(m_ulps);
            puts("m_output: "); mpfr_dump(m_output);
            puts("m_input_base: "); mpfr_dump(m_input_base);
            puts("m_input_exponent: "); mpfr_dump(m_input_exponent);
            puts("m_approx_error: "); mpfr_dump(m_approx_error);
            puts("m_ulps: "); mpfr_dump(m_ulps);
            puts("m_approx_error_ulps: "); mpfr_dump(m_approx_error_ulps);
            printf("i=%d  rounding_error=%g  input_base=%.17e  input_exponent=%.17e  approx_error_ulps=%g\n"
                   "nearest_fp64 =       %.17e\n"
                   "output_fp64_exact =  %.17e\n"
                   "output_fp64_approx = %.17e\n\n",
                   i, rounding_error, input_base_fp64, input_exponent_fp64, approx_error_ulps,
                   nearest_fp64, output_fp64_exact, output_fp64_approx);
        }

        approx_error_ulps_a[i] = fabs(approx_error_ulps);
        if (fabs(approx_error_ulps) > max_error_ulps) {
            max_error_ulps = fabs(approx_error_ulps);
        }
    }
    int unequal_count = N - equal_count;
    double unequal_ratio = unequal_count / (double) N;
    printf("%7d / %7d (%g%%) were bitwise unequal\n", unequal_count, N, 100 * unequal_ratio);
    printf("Max approximation error: %g ULPS\n", max_error_ulps);
    printf("Skipped ULPS computation for %d values where both approx and exact were rounded to zero\n", both_zero_count);
    printf("Skipped check for %d values where the exact answer was subnormal\n", subnormal_answer_count);

#if 0
    qsort(approx_error_ulps_a, N, sizeof(double), &compare_double);

    double error_percentile_25 = approx_error_ulps_a[(N*25L)/100];
    double error_percentile_50 = approx_error_ulps_a[(N*50L)/100];
    double error_percentile_75 = approx_error_ulps_a[(N*75L)/100];
    double error_percentile_99 = approx_error_ulps_a[(N*99L)/100];
    printf("25th percentile: %g ULPS\n", error_percentile_25);
    printf("50th percentile: %g ULPS\n", error_percentile_50);
    printf("75th percentile: %g ULPS\n", error_percentile_75);
    printf("99th percentile: %g ULPS\n", error_percentile_99);
    printf("Max: %g ULPS\n", approx_error_ulps_a[N-1]);

    for (int i = N-10; i < N; i++) {
        printf("approx_error_ulps_a[%d]: %g ULPS\n", i, approx_error_ulps_a[i]);
    }
#endif


    free(approx_error_ulps_a);

    mpfr_clear(m_input_base);
    mpfr_clear(m_input_exponent);
    mpfr_clear(m_output);
    mpfr_clear(m_rounding_error);
    mpfr_clear(m_approx_error);
    mpfr_clear(m_approx_error_ulps);
    mpfr_clear(m_ulps);
    mpfr_clear(m_nearest);
    mpfr_free_cache();
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

    snprintf(filename, filename_size, "%s_%s.bin", basename, "out_exp");
    load_binary_file_into_array(filename, exp_output, N);

    snprintf(filename, filename_size, "%s_%s.bin", basename, "out_expm1");
    load_binary_file_into_array(filename, expm1_output, N);

    snprintf(filename, filename_size, "%s_%s.bin", basename, "out_log");
    load_binary_file_into_array(filename, log_output, N);

    snprintf(filename, filename_size, "%s_%s.bin", basename, "out_pow");
    load_binary_file_into_array(filename, pow_output, N);

    // check output accuracy
    check_univariate_f(exp_input, exp_output, N, mpfr_exp, "exp");
    check_univariate_f(expm1_input, expm1_output, N, mpfr_expm1, "expm1");
    check_univariate_f(log_input, log_output, N, mpfr_log, "log");
    check_pow(pow_input_b, pow_input_e, pow_output, N);


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
