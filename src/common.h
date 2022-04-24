#ifndef _COMMON_H
#define _COMMON_H

#include "cellmodel.h"

#ifdef __cplusplus
extern "C" {
#endif

int dir_exists(char *dir);
void create_dir_if_missing(char *dir);
void memcpy_state(cellmodel_float_t *__restrict states_dst,
                  cellmodel_float_t *__restrict states_src, const struct cellmodel *model,
                  const long padded_num_cells, int state_index);
void memcpy_states(cellmodel_float_t *__restrict states_dst,
                   cellmodel_float_t *__restrict states_src, const struct cellmodel *model,
                   const long padded_num_cells);

#ifdef __cplusplus
}
#endif

#endif
