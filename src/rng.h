#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

uint64_t rng_draw_seed();
void add_random_noise(double *a, int N, uint64_t rng_seed, double amplitude);

#ifdef __cplusplus
}
#endif
