#ifndef _BENCH_H
#define _BENCH_H

#ifdef __cplusplus
extern "C" {
#endif

#include "bench_util.h"
#include "cellmodel.h"

void bench(const struct cellmodel *model, struct options *opt);

#ifdef __cplusplus
} // extern "C"
#endif

#endif
