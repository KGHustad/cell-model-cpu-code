#include "cellmodel.h"

void compare_naive_to_vanilla(const struct cellmodel *model_naive,
                              const struct cellmodel *model_vanilla, long num_cells, int num_steps,
                              scheme_type scheme, int show_rrms);
