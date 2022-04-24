extern "C" {
#include "bench_util.h"
}

#include "lut/LUT.hpp"

#include "cellmodel_lut.hpp"


void bench_lut(const struct cellmodel_lut *model, struct options *opt, default_LUT_type &lut_V,
               default_LUT_type &lut_Ca);
