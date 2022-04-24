#ifndef _TRACE_H
#define _TRACE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "cellmodel.h"

void dump_trace(cellmodel_float_t *t_trace, cellmodel_float_t *V_trace, int length, double solve_dt,
                scheme_type scheme, char *basename);

#ifdef __cplusplus
}
#endif

#endif
