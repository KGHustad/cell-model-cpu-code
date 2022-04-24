#ifndef _SOLVE_H
#define _SOLVE_H

#include "cellmodel.h"

void solve_single(const struct cellmodel *model, double dt, double T_end, int store_period,
                  scheme_type scheme, cellmodel_float_t **V_trace_ptr,
                  cellmodel_float_t **t_trace_ptr, int *trace_length);

#endif
