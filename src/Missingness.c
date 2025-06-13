
// File: Missingness.c
// Date: 13 May 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Mask out genotypes in replicates.

#include "Missingness.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "wavelib/wauxlib.h"
#include <time.h>
#include <math.h>

Mask_t* init_mask(int numSamples, int numRecords) {
    Mask_t* mask = (Mask_t*) calloc(1, sizeof(Mask_t));
    mask -> numSamples = numSamples;
    mask -> numRecords = numRecords;
    mask -> missing = (int**) calloc(numRecords, sizeof(int*));
    for (int i = 0; i < numRecords; i++)
        mask -> missing[i] = (int*) calloc(numSamples, sizeof(int));
    return mask;
}

// Fisher-Yates shuffle algorithm for an integer array.
void shuffle_real_array(gsl_rng* r, int* array, int n) {
    int temp = 0, j = 0;
    for (int i = n - 1; i > 0; i--) {
        j = gsl_rng_uniform_int(r, i + 1);
        temp = array[j];
        array[j] = array[i];
        array[i] = temp;
    }
}

MissingDistribution_t* init_missing_distribution(Replicate_t* replicate) {
    if (replicate == NULL)
        return NULL;
    if (replicate -> numRecords == 0 || replicate -> numSamples == 0)
        return NULL;

    MissingDistribution_t* dis = calloc(1, sizeof(MissingDistribution_t));
    
    dis -> proportions = (double*) calloc(replicate -> numRecords, sizeof(double));
    int curRecord = 0;
    // Calculate the proportion of missing sampels at each record.
    for (Record_t* temp = replicate -> headRecord; temp != NULL; temp = temp -> nextRecord) {
        int numMissing = 0;
        for (int i = 0; i < replicate -> numSamples; i++)
            if (temp -> genotypes[i].left == MISSING && temp -> genotypes[i].right == MISSING)
                numMissing++;
        dis -> proportions[curRecord] = numMissing /  (double) replicate -> numSamples;
        curRecord++;
    }
    dis -> numRecords = curRecord;
    
    return dis;
}

Mask_t* create_missing_mask(MissingDistribution_t* dis, int numSamples, int numRecords) {
    if (dis == NULL || numSamples == 0 || numRecords == 0)
        return NULL;

    // Create an empty mask.
    Mask_t* mask = init_mask(numSamples, numRecords);

    double* proportions = calloc(numRecords, sizeof(double));

    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));

    int* permu = (int*) calloc(numSamples, sizeof(int));
    for (int i = 0; i < numSamples; i++)
        permu[i] = i;

    if (numRecords <= dis -> numRecords) {
        
    } else {

    }

    for (int i = 0; i < numRecords; i++) {
        shuffle_real_array(r, permu, numSamples);
        int numMissing = (int) numSamples * proportions[i];
        for (int j = 0; j < numMissing; j++)
            mask -> missing[i][permu[j]] = MISSING;
    }

    free(proportions);
    gsl_rng_free(r);
    free(permu);

    return mask;
}

Mask_t* create_random_mask(MissingDistribution_t* dis, int numSamples, int numRecords, double mean, double stder) {

    // Create our mask.
    Mask_t* mask = init_mask(numSamples, numRecords);

    // Convert mean and stderr to alpha and beta that define a beta distribution.
    double alpha = -1;
    double beta = -1;
    if (dis == NULL) {
        alpha = (mean * mean * (1 - mean)) / (stder * stder) - mean;
        beta = (alpha / mean) * (1 - mean);
    }

    // Create our random number generator.
    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));

    // Create an array from 1 to numSamples.
    int* permu = (int*) calloc(numSamples, sizeof(int));
    for (int i = 0; i < numSamples; i++)
        permu[i] = i;

    // For each record in the partition.
    for (int i = 0; i < numRecords; i++) {
        // Create random permutation of samples.
        shuffle_real_array(r, permu, numSamples);
        // Draw a random proportion of missing samples from a the supplied distribution or beta distribution.
        int numMissing = 0; 
        if (dis == NULL)
            numMissing = (int) numSamples * gsl_ran_beta(r, alpha, beta);
        else
            numMissing = (int) numSamples * dis -> proportions[(int) (dis -> numRecords * gsl_rng_uniform(r))];

        // Set missing for samples at that record.
        for (int j = 0; j < numMissing; j++)
            mask -> missing[i][permu[j]] = MISSING;
    }
    
    gsl_rng_free(r);
    free(permu);

    return mask;
}

void apply_fill(Replicate_t* replicate, Mask_t* mask, int fill) {
    Record_t* prevMisRecord = NULL;
    Record_t* temp = NULL;
    int prevMisGenoPos = -1;
    
    // For each sample.
    for (int i = 0; i < mask -> numSamples; i++) { 
        temp = replicate -> headRecord;
        prevMisRecord = NULL;
        // Iterate through the records.
        for (int j = 0; j < mask -> numRecords; j++) {
            if (mask -> missing[j][i] == MISSING) {
                // If the previous missing genotype for the sample is <= fill, set genotypes inbetween the
                //  two sites to missing.
                if (prevMisRecord != NULL && (temp -> position - prevMisGenoPos) <= fill) {
                    for (Record_t* k = prevMisRecord -> nextRecord; k != temp; k = k -> nextRecord) {
                        k -> genotypes[i].isPhased = false;
                        k -> genotypes[i].left = MISSING;
                        k -> genotypes[i].right = MISSING;
                    }
                }
                prevMisGenoPos = temp -> position;
                prevMisRecord = temp;
            }
            temp = temp -> nextRecord;
        }
    }
}

void destroy_mask(Mask_t* mask) {
    if (mask == NULL)
        return;
    if (mask -> missing != NULL) {
        for (int i = 0; i < mask -> numRecords; i++)
            free(mask -> missing[i]);
        free(mask -> missing);
    }
    free(mask);
}

void destroy_missing_distribution(MissingDistribution_t* dis) {
    if (dis == NULL)
        return;
    if (dis -> proportions != NULL)
        free(dis -> proportions);
    free(dis);
}