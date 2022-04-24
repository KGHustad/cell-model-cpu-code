#ifndef _SOLVE_LUT_HPP
#define _SOLVE_LUT_HPP

#include "cellmodel_lut.hpp"
//#include "lut/LUT.hpp"

void solve_single_lut(const struct cellmodel_lut *model, double dt, double T_end, int store_period,
                      scheme_type scheme, cellmodel_float_t **V_trace_ptr,
                      cellmodel_float_t **t_trace_ptr, int *trace_length, default_LUT_type &lut_V,
                      default_LUT_type &lut_Ca);

#endif
