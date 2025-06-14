
// File: Interface.h
// Date: 13 May 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Set user parameters from CLI.

#ifndef INTERFACE_H
#define INTERFACE_H

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#include <stdbool.h>

// Hold CLI parameters.
typedef struct {
    // If records should be unphased.
    bool unphase;
    // If biallelic records should be unpolarized.
    bool unpolarize;
    // If records should be made into pseudohaps.
    bool pseudohap;
    // The out file basename.
    char* outFile;
    // The mask file name.
    char* maskFile;
    // The string defining the VCF or mean and std. dev. for beta distribution.
    char* betaMissing;
    // The VCF file for to calculate the distribution of missing genotypes.
    char* randomMissing;
    // The mean and std err for random missing genotypes.
    double meanMissing;
    double stdMissing;
    // The segment length in base pairs for ms replicates.
    int length;
    // The EGGS command the user ran.
    char* command;
} EggsConfig_t;

// Accepts the command line arguments and sets configuration parameters.
// Accepts:
//  int argc, char *argv[] -> The command line arguments.
// Returns: EggsConfig_t*, the configuration. NULL, if invalid parameter given or CLI parsing error.
EggsConfig_t* init_eggs_configuration(int argc, char *argv[]);

// Frees memory allocated to configuration.
void destroy_eggs_configuration(EggsConfig_t* eggsConfig);

// Prints help menu.
void print_help();

#endif