
// File: Interface.c
// Date: 13 May 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Set user parameters from CLI.

#include "Interface.h"
#include <stdio.h>
#include "ketopt.h"
#include <unistd.h>
#include "kstring.h"

// Our help menu.
void print_help() {
    fprintf(stderr, "\n");
    fprintf(stderr, "EGGS v1.0\n");
    fprintf(stderr, "---------\n\n");
    fprintf(stderr, "Written by T. Quinn Smith\n");
    fprintf(stderr, "Principal Investigator: Zachary A. Szpiech\n");
    fprintf(stderr, "The Pennsylvania State University\n\n");
    fprintf(stderr, "Usage: eggs [OPTIONS]\n\n");
    fprintf(stderr, "Reads from stdin and write to stdout unless -o is provided.\n");
    fprintf(stderr, "OPTIONS:\n");
    fprintf(stderr, "    -h,--help                       Print help.\n");
    fprintf(stderr, "    -u,--unphase                    Left and right genotypes are swapped with a probability of 0.5\n");
    fprintf(stderr, "    -p,--unpolarize                 Biallelic site alleles swapped with a probability of 0.5\n");
    fprintf(stderr, "    -s,--pseudohap                  Pseudohaploidize all samples. Automatically removes phase.\n");
    fprintf(stderr, "    -o,--out        STR             Basename to use for output files instead of stdout.\n");
    fprintf(stderr, "    -m,--mask       VCF             Filename of VCF to use as mask for missing genotypes.\n");
    fprintf(stderr, "    -b,--beta       VCF/STR         Calculate mu/sigma of missingness per site from VCF or supply as\n");
    fprintf(stderr, "                                        values as \"mu,sigma\". Defines beta dsitribution for missingness.\n");
    fprintf(stderr, "    -r,--random     VCF             Calculates proportion of missing samples per site from VCF and uses that\n");
    fprintf(stderr, "                                        distribution to randomly introduce missing genotypes.\n");
    fprintf(stderr, "    -d,--deamin     STR             Two comma-seperated proportions \"prop1,prop2\" where prop1 is the probability\n");
    fprintf(stderr, "                                        the site is a transition and prop2 is the probability of deamination.\n");
    fprintf(stderr, "    -l,--length     INT             Only used with ms-style input. Sets length of segments in base-pairs.\n");
    fprintf(stderr, "                                        Default 1,000,000 base-pairs.\n");
    fprintf(stderr, "    -a,--hap                        Only used with ms-style input. Each diploid only has one lineage.\n");
    fprintf(stderr, "\n");
}

// Our options.
static ko_longopt_t long_options[] = {
    {"help",            ko_no_argument,         'h'},
    {"unphase",         ko_no_argument,         'u'},
    {"unpolarize",      ko_no_argument,         'p'},
    {"pseudohap",       ko_no_argument,         's'},
    {"out",             ko_required_argument,   'o'},
    {"mask",            ko_required_argument,   'm'},
    {"beta",            ko_required_argument,   'b'},
    {"random",          ko_required_argument,   'r'},
    {"length",          ko_required_argument,   'l'},
    {"hap",             ko_no_argument,         'a'},
    {"deamin",          ko_required_argument,   'b'},
    {0, 0, 0}
};

// Accepts parsed parameters from user and ensures they are valid.
// Returns: 0, if valid. -1, if invalid.
int check_configuration(EggsConfig_t* eggsConfig) {
    // If length was set by the user, it must be an integer >= 1000.
    if (eggsConfig -> length != -1 && eggsConfig -> length < 1000) {
        fprintf(stderr, "-l must be given an integer >= 1000. Exiting!\n");
        return -1;
    }
    // Cannot use mask, beta, or random genotypes together.
    int numSet = (int) (eggsConfig -> randomMissing != NULL) + (int) (eggsConfig -> maskFile != NULL) + (int) (eggsConfig -> betaMissing != NULL);
    if (numSet > 1) {
        fprintf(stderr, "Cannot use -m, -b, and -r options together. Exiting!\n");
        return -1;
    }
    // If mask file given, make sure it exists.
    if (eggsConfig -> maskFile != NULL && access(eggsConfig -> maskFile, F_OK) != 0) {
        fprintf(stderr, "-m %s does not exist. Exiting!\n", eggsConfig -> maskFile);
        return -1;
    }
    // If beta file given, make sure it exists.
    if (eggsConfig -> randomMissing != NULL && access(eggsConfig -> randomMissing, F_OK) != 0) {
        fprintf(stderr, "-r %s does not exist. Exiting!\n", eggsConfig -> randomMissing);
        return -1;
    }
    // If beta was given and VCF file does not exists, then parser values directly.
    if (eggsConfig -> betaMissing != NULL && access(eggsConfig -> betaMissing, F_OK) != 0) {
        char* meanstd = strdup(eggsConfig -> betaMissing);
        char* next  = NULL;
        eggsConfig -> meanMissing = strtod(meanstd, &next);
        // Ensure mean and stder were given seperated by a comma.
        if (meanstd == next || next[0] != ',') {
            fprintf(stderr, "-r must be given two positive real numbers seperated by a comma. Exiting!\n");
            free(meanstd);
            return -1;
        }
        eggsConfig -> stdMissing = strtod(next + 1, (char**) NULL);
        // Make sure mean and stder are valid for a beta distribution.
        if (eggsConfig -> meanMissing >= 1 || eggsConfig -> meanMissing <= 0 || eggsConfig -> stdMissing <= 0 || eggsConfig -> stdMissing * eggsConfig -> stdMissing >= eggsConfig -> meanMissing * (1 - eggsConfig -> meanMissing)) {
            fprintf(stderr, "-r must satisfy parameters for a beta distribution. Exiting!\n");
            free(meanstd);
            return -1;
        }
        free(meanstd);
    }
    // If deamination parameters where given.
    if (eggsConfig -> deamin != NULL) {
        char* deamin = strdup(eggsConfig -> deamin);
        char* next  = NULL;
        eggsConfig -> probTransition = strtod(deamin, &next);
        // Ensure both proportions were given seperated by a comma.
        if (deamin == next || next[0] != ',') {
            fprintf(stderr, "-d must be given two positive real numbers seperated by a comma. Exiting!\n");
            free(deamin);
            return -1;
        }
        eggsConfig -> probDeamination = strtod(next + 1, (char**) NULL);
        if (eggsConfig -> probDeamination >= 1 || eggsConfig -> probDeamination <= 0 || eggsConfig -> probTransition >= 1 || eggsConfig -> probTransition <= 0) {
            fprintf(stderr, "-d must be given two real numbers in (0, 1). Exiting!\n");
            free(deamin);
            return -1;
        }
        free(deamin);
    }
    return 0;
}

EggsConfig_t* init_eggs_configuration(int argc, char *argv[]) {

    const char *opt_str = "haupso:m:r:l:b:d:";
    ketopt_t options = KETOPT_INIT;
    int c;

    // Make sure we can parse the options. If -h was given, print help menu.
    while ((c = ketopt(&options, argc, argv, 1, opt_str, long_options)) >= 0) {
        switch (c) {
            case ':': fprintf(stderr, "Error! Option %s is missing an argument! Exiting ...\n", argv[options.i - 1]); return NULL;
            case '?': fprintf(stderr, "Error! \"%s\" is unknown! Exiting ...\n", argv[options.i - 1]); return NULL;
            case 'h': print_help(); return NULL;
        }
	}

    // Set configuration defaults.
    EggsConfig_t* eggsConfig = (EggsConfig_t*) calloc(1, sizeof(EggsConfig_t));
    eggsConfig -> unphase = false;
    eggsConfig -> unpolarize = false;
    eggsConfig -> pseudohap = false;
    eggsConfig -> outFile = NULL;
    eggsConfig -> maskFile = NULL;
    eggsConfig -> betaMissing = NULL;
    eggsConfig -> randomMissing = NULL;
    eggsConfig -> meanMissing = -1;
    eggsConfig -> stdMissing = -1;
    eggsConfig -> length = 1000000;
    eggsConfig -> hap = false;
    eggsConfig -> deamin = NULL;
    eggsConfig -> probDeamination = 0;
    eggsConfig -> probTransition = 0;
    eggsConfig -> command = NULL;

    // Get parameters from user.
    options = KETOPT_INIT;
    while ((c = ketopt(&options, argc, argv, 1, opt_str, long_options)) >= 0) {
        switch (c) {
            case 'u': eggsConfig -> unphase = true; break;
            case 'p': eggsConfig -> unpolarize = true; break;
            case 's': eggsConfig -> pseudohap = true; break;
            case 'o': eggsConfig -> outFile = strdup(options.arg); break;
            case 'm': eggsConfig -> maskFile = strdup(options.arg); break;
            case 'r': eggsConfig -> randomMissing = strdup(options.arg); break;
            case 'b': eggsConfig -> betaMissing = strdup(options.arg); break;
            case 'l': eggsConfig -> length = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 'a': eggsConfig -> hap = true;
            case 'd': eggsConfig -> deamin = strdup(options.arg); break;
        }
	}

    // If a parameter was invalid, we return null.
    if (check_configuration(eggsConfig) != 0) {
        destroy_eggs_configuration(eggsConfig);
        return NULL;
    }

    // Copy user EGGS command.
    kstring_t* cmd = (kstring_t*) calloc(1, sizeof(kstring_t));
    for (int i = 0; i < argc; i++) 
        ksprintf(cmd, "%s ", argv[i]);
    eggsConfig -> command = cmd -> s;
    free(cmd);

    return eggsConfig;
}

void destroy_eggs_configuration(EggsConfig_t* eggsConfig) {
    if (eggsConfig == NULL)
        return;
    if (eggsConfig -> outFile != NULL)
        free(eggsConfig -> outFile);
    if (eggsConfig -> maskFile != NULL)
        free(eggsConfig -> maskFile);
    if (eggsConfig -> randomMissing != NULL)
        free(eggsConfig -> randomMissing);
    if (eggsConfig -> betaMissing != NULL)
        free(eggsConfig -> betaMissing);
    if (eggsConfig -> command != NULL)
        free(eggsConfig -> command);
    if (eggsConfig -> deamin != NULL)
        free(eggsConfig -> deamin);
    free(eggsConfig);
}