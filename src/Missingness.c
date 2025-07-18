
// File: Missingness.c
// Date: 13 May 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Mask out genotypes in replicates.

#include "Missingness.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include <time.h>
#include <math.h>

double uniform(gsl_rng* r, double min, double max) {
    return min + (max - min) * (double) gsl_rng_uniform(r);
}

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

MissingDistribution_t* init_missing_distribution(Replicate_t* replicate, InputStream_t* inputStream) {
    if (replicate == NULL)
        return NULL;
    if (inputStream == NULL)
        return NULL;
    if (replicate -> numSamples == 0)
        return NULL;

    MissingDistribution_t* dis = calloc(1, sizeof(MissingDistribution_t));

    typedef struct Proportions{
        double prop;
        struct Proportions* next;
    } Proportions_t;

    Proportions_t* head = NULL;
    Proportions_t* tail = NULL;
    Proportions_t* temp = NULL;
    
    Record_t* record = (Record_t*) calloc(1, sizeof(Record_t));
    record -> genotypes = (Genotype_t*) calloc(replicate -> numSamples, sizeof(Genotype_t));
    record -> numSamples = replicate -> numSamples;
    
    int numRecords = 0;
    while (get_next_vcf_record(record, inputStream)) {
        int numMissing = 0;
        for (int i = 0; i < replicate -> numSamples; i++)
            if (record -> genotypes[i].left == MISSING && record -> genotypes[i].right == MISSING)
                numMissing++;
        temp = calloc(1, sizeof(Proportions_t*));
        if (numRecord == 0) {
            temp -> prop = numMissing / (double) replicate -> numSamples;
            head = temp;
            tail = temp;
        } else {
            tail -> next = temp;
            tail = temp;
        }
        numRecords++;
    }

    dis -> proportions = (double*) calloc(numRecords, sizeof(double));
    dis -> numRecords = numRecords;
    temp = head;
    for (int i = 0; i < numRecords; i++) {
        dis -> proportions[i] = temp -> prop;
        temp = head;
        head = head -> next;
        free(temp);
    }

    destroy_record(record);
    
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

    // If we are shrinking the dispersal
    if (numRecords <= dis -> numRecords) {
        int chunkSize = (int) ceil(dis -> numRecords / (double) numRecords);
        int nextChunk = 0;
        // Randomly choose a proportion of missingess from the chunk of sites.
        for (int i = 0; i < numRecords - 1; i++) {
            proportions[i] = dis -> proportions[nextChunk + gsl_rng_uniform_int(r, chunkSize)];
            nextChunk += chunkSize;
        }
        proportions[numRecords - 1] = dis -> proportions[nextChunk + gsl_rng_uniform_int(r, dis -> numRecords - nextChunk)];
    // If we are stretching the dispersal
    } else {
        int chunkSize = (int) ceil(numRecords / (double) dis -> numRecords);
        int nextChunk = 0;
        // For each chunk, take two adjacent sites and randomly get a proportion between the two bounds.
        for (int i = 0; i < dis -> numRecords - 1; i++) {
            for (int j = 0; j < chunkSize; j++)
                 proportions[nextChunk + j] = uniform(r, fmin(dis -> proportions[i], dis -> proportions[i + 1]), fmax(dis -> proportions[i], dis -> proportions[i + 1]));
            nextChunk += chunkSize;
        }
        // We model the last chunk after the last two sites.
        for (int j = nextChunk; j < numRecords; j++)
            proportions[j] = uniform(r, fmin(dis -> proportions[dis -> numRecords - 2], dis -> proportions[dis -> numRecords - 1]), fmax(dis -> proportions[dis -> numRecords - 2], dis -> proportions[dis -> numRecords - 1]));
    }

    // Randomly introduce missingness to the samples at each site.
    for (int i = 0; i < numRecords; i++) {
        if (proportions[i] == 0)
            continue;
        else if (proportions[i] == 1)
            for (int j = 0; j < numSamples; j++)
                mask -> missing[i][permu[j]] = MISSING;
        else {
            shuffle_real_array(r, permu, numSamples);
            int numMissing = (int) numSamples * proportions[i];
            for (int j = 0; j < numMissing; j++)
                mask -> missing[i][permu[j]] = MISSING;
        }
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