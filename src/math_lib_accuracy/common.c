#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "common.h"

int load_init_file(char *basename, int *N_ptr)
{
    int filename_size = 1024;
    char filename[filename_size];
    snprintf(filename, filename_size, "%s_init.txt", basename);

    FILE *f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "Error: Could not open init file at '%s'\n", filename);
        exit(EXIT_FAILURE);
    }


    int N;
    int read = fscanf(f, "%d", &N);
    if (read != 1) {
        fprintf(stderr, "Error: Expected integer at the start of the file '%s'\n", filename);
        exit(EXIT_FAILURE);
    }

    assert(N > 0);

    fclose(f);

    *N_ptr = N;
    return 0;
}

int load_binary_file_into_array(char *filename, double *b, int N)
{
    assert(b != NULL);

    FILE *f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "Error: Could not open binary file at '%s' for reading\n", filename);
        exit(EXIT_FAILURE);
    }

    size_t elements_read = fread(b, sizeof(double), N, f);
    assert(elements_read == N);

    fclose(f);

    return 0;
}
