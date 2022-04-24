#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "schemes.h"
#include "trace.h"

void dump_trace(cellmodel_float_t *t_trace, cellmodel_float_t *V_trace, int length, double solve_dt,
                scheme_type scheme, char *basename)
{
    const size_t outfilename_size = 1024;
    char outfilename[outfilename_size];
    assert(strlen(basename) < 1000);

    FILE *f;

    snprintf(outfilename, outfilename_size, "%s_V_trace.bin", basename);
    f = fopen(outfilename, "w");
    assert(f != NULL);
    fwrite(V_trace, sizeof(cellmodel_float_t), length, f);
    fclose(f);

    snprintf(outfilename, outfilename_size, "%s_t_trace.bin", basename);
    f = fopen(outfilename, "w");
    assert(f != NULL);
    fwrite(t_trace, sizeof(cellmodel_float_t), length, f);
    fclose(f);

    snprintf(outfilename, outfilename_size, "%s_info.txt", basename);
    f = fopen(outfilename, "w");
    assert(f != NULL);
    {
        size_t scheme_str_size = 100;
        char scheme_str[scheme_str_size];
        get_scheme_str(scheme_str, scheme_str_size, scheme);
        fprintf(f, "scheme=%s\n", scheme_str);
    }
    fprintf(f, "solve_dt=%g\n", solve_dt);
    fclose(f);
}
