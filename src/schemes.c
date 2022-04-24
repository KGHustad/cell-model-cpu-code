#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "schemes.h"

void get_scheme_str(char *scheme_str, size_t str_size, scheme_type scheme)
{
    switch (scheme) {
    case SCHEME_FE:
        snprintf(scheme_str, str_size, "FE");
        break;
    case SCHEME_RL:
        snprintf(scheme_str, str_size, "RL");
        break;
    case SCHEME_GRL1:
        snprintf(scheme_str, str_size, "GRL1");
        break;
    default:
        printf("ERROR: Scheme %d has no known name\n", scheme);
        exit(EXIT_FAILURE);
    }
}
