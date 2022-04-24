#include "cellmodel.h"
#include "cellmodel_lut.hpp"

void compare_naive_to_lut(const struct cellmodel *model_naive,
                          const struct cellmodel_lut *model_lut, long num_cells, int num_steps,
                          double V_step, double Ca_step, scheme_type scheme);
