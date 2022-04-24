#include <random>

#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include "rng.h"

extern "C" uint64_t rng_draw_seed()
{
    FILE *f = fopen("/dev/urandom", "r");
    uint64_t rng_seed;
    int items_read = fread(&rng_seed, sizeof(rng_seed), 1, f);
    assert(items_read == 1);
    fclose(f);
    return rng_seed;
}

extern "C" void add_random_noise(double *a, int N, uint64_t rng_seed, double amplitude)
{
    std::mt19937_64 rng;
    rng.seed(rng_seed);

    std::uniform_real_distribution<double> dist(-amplitude, amplitude);
    for (int i = 0; i < N; i++) {
        a[i] += dist(rng);
    }
}
