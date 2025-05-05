
#include "Interface.h"

#include <stdio.h>
#include <ketopt.h>
#include <unistd.h>

void print_help() {
    printf("\n");
    printf("EGGS v1.0 May 2025\n");
    printf("----------------------\n\n");
    printf("Written by T. Quinn Smith\n");
    printf("Principal Investigator: Zachary A. Szpiech\n");
    printf("The Pennsylvania State University\n\n");
    printf("Usage: eggs [options]\n\n");
    printf("Options:\n");
    printf("    -h,--help                       Print help.\n");
    printf("    -u,--unphase                    Left and right genotypes are swapped with a probability of 0.5\n");
    printf("    -p,--unpolarize                 Biallelic site alleles swapped with a probability of 0.5\n");
    printf("    -s,--pseudohap                  Pseudohaplotize all samples.\n");
    printf("    -o,--out        STR             Basename to use for output files instead of stdout.\n");
    printf("    -m,--mask       STR             Filename of VCF to use as mask for missing genotypes.\n");
    printf("    -f,--fill       INT             Used with -m. If distance (in base-paris) between missing\n");
    printf("                                        sites is <= INT, then sample's genotypes between are set to missing.\n");
    printf("    -r,--random     DOUBLE,DOUBLE   The mean and standard error used to introduce missing genotypes.\n");
    printf("    -l,--length     INT             Only used with ms-style input. Sets length of segments in base-pairs.\n");
    printf("                                        Default 1,000,000 base-pairs.\n");
    printf("\n");
}

static ko_longopt_t long_options[] = {
    {"help",            ko_no_argument,         'h'},
    {"unphase",         ko_no_argument,         'u'},
    {"unpolarize",      ko_no_argument,         'p'},
    {"pseudohap",       ko_no_argument,         's'},
    {"out",             ko_required_argument,   'o'},
    {"mask",            ko_required_argument,   'm'},
    {"fill",            ko_required_argument,   'f'},
    {"random",          ko_required_argument,   'r'},
    {"length",          ko_required_argument,   'l'},
    {0, 0, 0}
};

int check_configuration(EggsConfig_t* eggsConfig) {
    if (eggsConfig -> fill != -1 && eggsConfig -> fill <= 0) {
        fprintf(stderr, "-f must be given an integer > 0. Exiting!\n");
        destroy_eggs_configuration(eggsConfig);
        return -1;
    }
    if (eggsConfig -> length != -1 && eggsConfig -> length < 1000) {
        fprintf(stderr, "-l must be given an integer >= 1000. Exiting!\n");
        destroy_eggs_configuration(eggsConfig);
        return -1;
    }
    if (eggsConfig -> maskFile != NULL && access(eggsConfig -> maskFile, F_OK) != 0) {
        fprintf(stderr, "-m %s does not exist. Exiting!\n", eggsConfig -> maskFile);
        destroy_eggs_configuration(eggsConfig);
        return -1;
    }
    if (eggsConfig -> randomMissing != NULL) {
        char* meanstd = strdup(eggsConfig -> randomMissing);
        char* next  = NULL;
        eggsConfig -> meanMissing = strtod(meanstd, &next);
        if (meanstd == next) {
            fprintf(stderr, "-r must be given two positive real numbers seperated by a comma. Exiting!\n");
            destroy_eggs_configuration(eggsConfig);
            free(meanstd);
            return -1;
        }
        eggsConfig -> stdMissing = strtod(next + 1, (char**) NULL);
        if (eggsConfig -> meanMissing < 0 || eggsConfig -> stdMissing <= 0) {
            fprintf(stderr, "-r must be given a mean >= 0 and a stderr > 0 seperated by a comma. Exiting!\n");
            free(meanstd);
            destroy_eggs_configuration(eggsConfig);
            return -1;
        }
        free(meanstd);
    }
    return 0;
}

EggsConfig_t* init_eggs_configuration(int argc, char *argv[]) {

    const char *opt_str = "hupso:m:f:r:l:";
    ketopt_t options = KETOPT_INIT;
    int c;

    while ((c = ketopt(&options, argc, argv, 1, opt_str, long_options)) >= 0) {
        switch (c) {
            case ':': fprintf(stderr, "Error! Option %s is missing an argument! Exiting ...\n", argv[options.i - 1]); return NULL;
            case '?': fprintf(stderr, "Error! \"%s\" is unknown! Exiting ...\n", argv[options.i - 1]); return NULL;
            case 'h': print_help(); return NULL;
        }
	}

    EggsConfig_t* eggsConfig = calloc(1, sizeof(EggsConfig_t));
    eggsConfig -> unphase = false;
    eggsConfig -> unpolarize = false;
    eggsConfig -> pseudohap = false;
    eggsConfig -> outFile = NULL;
    eggsConfig -> maskFile = NULL;
    eggsConfig -> fill = -1;
    eggsConfig -> meanMissing = -1;
    eggsConfig -> stdMissing = -1;
    eggsConfig -> length = 1000000;

    options = KETOPT_INIT;
    while ((c = ketopt(&options, argc, argv, 1, opt_str, long_options)) >= 0) {
        switch (c) {
            case 'u': eggsConfig -> unphase = true; break;
            case 'p': eggsConfig -> unpolarize = true; break;
            case 's': eggsConfig -> pseudohap = true; break;
            case 'o': eggsConfig -> outFile = strdup(options.arg); break;
            case 'm': eggsConfig -> maskFile = strdup(options.arg); break;
            case 'f': eggsConfig -> fill = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 'r': eggsConfig -> randomMissing = strdup(options.arg); break;
            case 'l': eggsConfig -> length = (int) strtol(options.arg, (char**) NULL, 10); break;
        }
	}

    if (check_configuration(eggsConfig) != 0)
        return NULL;

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
    free(eggsConfig);
}