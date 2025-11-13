
// File: Missingness.c
// Date: 13 May 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Mask out genotypes in replicates.

#include "Missingness.h"
#include "kvec.h"
#include <time.h>
#include <math.h>

double uniform(gsl_rng* r, double min, double max) {
    return min + (max - min) * (double) gsl_rng_uniform(r);
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

    kvec_t(double) proportions;
    kv_init(proportions);

    Record_t* record = (Record_t*) calloc(1, sizeof(Record_t));
    record -> genotypes = (Genotype_t*) calloc(replicate -> numSamples, sizeof(Genotype_t));
    record -> numSamples = replicate -> numSamples;
    
    int numRecords = 0;
    while (get_next_vcf_record(record, inputStream, false)) {
        int numMissing = 0;
        for (int i = 0; i < replicate -> numSamples; i++)
            if (record -> genotypes[i].left == MISSING && record -> genotypes[i].right == MISSING)
                numMissing++;
        kv_push(double, proportions, numMissing / (double) record -> numSamples);
        numRecords++;
    }

    dis -> proportions = (double*) calloc(numRecords, sizeof(double));
    dis -> numRecords = numRecords;
    for (int i = 0; i < numRecords; i++) 
        dis -> proportions[i] = kv_A(proportions, i);

    destroy_record(record);
    kv_destroy(proportions);

    // Create RNG.
    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));
    dis -> r = r;

    dis -> permu = NULL;
    
    return dis;
}


void create_random_mask(MissingDistribution_t* dis, int numSamples, int curRecord, int numRecords, double mean, double stder, int* mask) {

    // If the shuffle array has not been created, then allocate array.
    if (dis -> permu == NULL) {
        dis -> permu = calloc(numSamples, sizeof(int));
        dis -> sizeOfPermu = numSamples;
        for (int i = 0; i < numSamples; i++)
            dis -> permu[i] = i;
    } else if (dis -> sizeOfPermu != numSamples) {
        free(dis -> permu);
        dis -> permu = calloc(numSamples, sizeof(int));
        dis -> sizeOfPermu = numSamples;
        for (int i = 0; i < numSamples; i++)
            dis -> permu[i] = i;
    }

    // Compute our number of missing samples.
    int numMissing;
    // Use beta-distribution.
    if (mean != -1 && stder != -1) {
        double alpha = (mean * mean * (1 - mean)) / (stder * stder) - mean;
        double beta = (alpha / mean) * (1 - mean);
        numMissing = (int) (numSamples * gsl_ran_beta(dis -> r, alpha, beta));
    // Pick randomly.
    } else if (curRecord == 0) {
        numMissing = (int) (numSamples * dis -> proportions[(int) (dis -> numRecords * gsl_rng_uniform(dis -> r))]);
    // Use EGGS method.
    } else {
        if (numRecords < dis -> numRecords) {
            int chunkSize = (int) ceil(dis -> numRecords / (double) numRecords);
            int chunkStart = chunkSize * (curRecord / chunkSize);
            if (chunkStart + chunkSize > dis -> numRecords) 
                numMissing = (int) (numSamples * dis -> proportions[chunkStart + gsl_rng_uniform_int(dis -> r, dis -> numRecords - chunkSize)]);
            else 
                numMissing = (int) (numSamples * dis -> proportions[chunkStart + gsl_rng_uniform_int(dis -> r, chunkSize)]);
        } else if (numRecords == dis -> numRecords) {
            numMissing = (int) (numSamples * dis -> proportions[curRecord]);
        } else {
            int chunkSize = (int) ceil(numRecords / (double) dis -> numRecords);
            int chunkStart = curRecord / chunkSize;
            if (chunkStart == dis -> numRecords - 1) 
                numMissing = (int) (numSamples * uniform(dis -> r, fmin(dis -> proportions[dis -> numRecords - 2], dis -> proportions[dis -> numRecords - 1]), fmax(dis -> proportions[dis -> numRecords - 2], dis -> proportions[dis -> numRecords - 1])));
            else
                numMissing = (int) (numSamples * uniform(dis -> r, fmin(dis -> proportions[chunkStart], dis -> proportions[chunkStart + 1]), fmax(dis -> proportions[chunkStart], dis -> proportions[chunkStart + 1])));
        }
    }

    // Set numSamples to missing.
    shuffle_real_array(dis -> r, dis -> permu, numSamples);
    for (int i = 0; i < numSamples; i++)
        mask[i] = 0;
    for (int i = 0; i < numMissing; i++) 
        mask[dis -> permu[i]] = MISSING;

}

void destroy_missing_distribution(MissingDistribution_t* dis) {
    if (dis == NULL)
        return;
    if (dis -> proportions != NULL)
        free(dis -> proportions);
    if (dis -> r != NULL)
        gsl_rng_free(dis -> r);
    if (dis -> permu != NULL)
        free(dis -> permu);
    free(dis);
}