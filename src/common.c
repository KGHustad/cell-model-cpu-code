#include "common.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

int dir_exists(char *dir)
{
    struct stat buf;
    return stat(dir, &buf) == 0 && S_ISDIR(buf.st_mode);
}

void create_dir_if_missing(char *dir)
{
    if (dir_exists(dir)) {
        return;
    }

    int retval = mkdir(dir, 0700);
    if (retval) {
        printf("ERROR: Failed to create directory '%s'.\n", dir);
        exit(EXIT_FAILURE);
    }
}


void memcpy_state(cellmodel_float_t *__restrict states_dst,
                  cellmodel_float_t *__restrict states_src, const struct cellmodel *model,
                  const long padded_num_cells, int state_index)
{
    const int s = state_index;
#pragma omp parallel for
    for (long i = 0; i < padded_num_cells; i++) {
        states_dst[s * padded_num_cells + i] = states_src[s * padded_num_cells + i];
    }
}

void memcpy_states(cellmodel_float_t *__restrict states_dst,
                   cellmodel_float_t *__restrict states_src, const struct cellmodel *model,
                   const long padded_num_cells)
{
    const int num_states = model->num_states;

#pragma omp parallel
    for (int s = 0; s < num_states; s++) {
#pragma omp for
        for (long i = 0; i < padded_num_cells; i++) {
            states_dst[s * padded_num_cells + i] = states_src[s * padded_num_cells + i];
        }
    }
}
