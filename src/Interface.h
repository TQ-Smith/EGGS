#ifndef INTERFACE_H
#define INTERFACE_H

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#include <stdbool.h>

typedef struct {
    bool unphase;
    bool unpolarize;
    bool pseudohap;
    char* outFile;
    char* maskFile;
    int fill;
    char* randomMissing;
    double meanMissing;
    double stdMissing;
    int length;
} EggsConfig_t;

EggsConfig_t* init_eggs_configuration(int argc, char *argv[]);

void destroy_eggs_configuration(EggsConfig_t* eggsConfig);

void print_help();

#endif