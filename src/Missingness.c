
// File: Missingness.c
// Date: 13 May 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Mask out genotypes in replicates.

#include "Missingness.h"
#include <time.h>
#include <math.h>

double uniform(gsl_rng* r, double min, double max) {
    return min + (max - min) * (double) gsl_rng_uniform(r);
}

void shuffle_real_array(gsl_rng* r, int* array, int n) {
    int temp = 0, j = 0;
    for (int i = n - 1; i > 0; i--) {
        j = gsl_rng_uniform_int(r, i + 1);
        temp = array[j];
        array[j] = array[i];
        array[i] = temp;
    }
}

MissingMask_t* init_missing_mask(Replicate_t* replicate, InputStream_t* inputStream) {
    if (replicate == NULL)
        return NULL;
    if (inputStream == NULL)
        return NULL;
    if (replicate -> numSamples == 0)
        return NULL;
    
    MissingMask_t* mask = calloc(1, sizeof(MissingMask_t));
    mask -> numSamples = replicate -> numSamples;
    mask -> blockMissing = calloc(replicate -> numSamples, sizeof(double));
    kv_init(mask -> mask);

    Record_t* record = (Record_t*) calloc(1, sizeof(Record_t));
    record -> genotypes = (Genotype_t*) calloc(replicate -> numSamples, sizeof(Genotype_t));
    record -> numSamples = replicate -> numSamples;

    // Read in each record and set bit corresponding to sample if missing.
    while (get_next_vcf_record(record, inputStream, false, false)) {
        CompactBitset* cb = cb_create(replicate -> numSamples);
        for (int i = 0; i < replicate -> numSamples; i++)
            if (record -> genotypes[i].left == MISSING && record -> genotypes[i].right == MISSING)
                cb_set_bit(cb, i);
        kv_push(CompactBitset*, mask -> mask, cb);
        mask -> numRecords++;
    }
    destroy_record(record);

    // Create RNG.
    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));
    mask -> r = r;
    
    return mask;
}

void get_mask_for_next_site(MissingMask_t* mask, CompactBitset* cb, int numRecords, int numSamples, int site) {
    // Allocate memory.
    if (mask -> prevCorrespondingSample == NULL && mask -> correspondingSample == NULL) {
        mask -> prevCorrespondingSample = calloc(numSamples, sizeof(int));
        mask -> correspondingSample = calloc(numSamples, sizeof(int));
    }

    // Get block boundaries.
    int lower = (int) (site * (mask -> numRecords / (double) numRecords));
    int upper;
    if (site == numRecords - 1) 
        upper = numRecords - 1;
    else 
        upper = (int) ((site + 1) * (mask -> numRecords / (double) numRecords)) - 1;
    int size = (upper - lower) + 1;

    // Calculate proportion of missing sites for samples within the current block.
    for (int i = 0; i < mask -> numSamples; i++)
        mask -> blockMissing[i] = 0;
    for (int l = lower; l <= upper; l++) {
        for (int i = 0; i < mask -> numSamples; i++) {
            CompactBitset* cb = kv_A(mask -> mask, l);
            if (cb_get_bit(cb, i))
                mask -> blockMissing[i]++;
        }
    }
    for (int i = 0; i < mask -> numSamples; i++)
        mask -> blockMissing[i] /= size;
    
    // If it is the first site, we randomly choose.
    if (site == 0) {
        for (int i = 0; i < numSamples; i++) {
            // Choose a random sample.
            int randSample = gsl_rng_uniform_int(mask -> r, mask -> numSamples);
            // With that sample's probability, determine if it should be missing.
            if ((double) gsl_rng_uniform(mask -> r) < mask -> blockMissing[randSample]) {
                // Set bit and save the random sample.
                cb_set_bit(cb, i);
                mask -> correspondingSample[i] = randSample;
            } else {
                // Otherwise no sample was chosen.
                mask -> correspondingSample[i] = -1;
            }
        }
    } else {
        for (int i = 0; i < numSamples; i++) {
            // If the previous site was missing. We use -1 as a flag to denote the previous site was not missing.
            if (mask -> prevCorrespondingSample[i] != -1) {
                // Use the same sample's probability if the current block should be missing.
                if ((double) gsl_rng_uniform(mask -> r) < mask -> blockMissing[mask -> prevCorrespondingSample[i]]) {
                    // Set bit.
                    cb_set_bit(cb, i);
                    // Keep current sample.
                    mask -> correspondingSample[i] = mask -> prevCorrespondingSample[i];
                } else {
                    mask -> correspondingSample[i] = -1;
                }
            // If the previous site was not missing, pick random sample.
            } else {
                int randSample = gsl_rng_uniform_int(mask -> r, mask -> numSamples);
                if ((double) gsl_rng_uniform(mask -> r) < mask -> blockMissing[randSample]) {
                    cb_set_bit(cb, i);
                    mask -> correspondingSample[i] = randSample;
                } else {
                    mask -> correspondingSample[i] = -1;
                }
            }
        }
    }
    // Swap sample tracking.
    int* temp = mask -> prevCorrespondingSample;
    mask -> prevCorrespondingSample = mask -> correspondingSample;
    mask -> correspondingSample = temp;
}

void destroy_missing_mask(MissingMask_t* mask) {
    if (mask == NULL)
        return;
    if (mask -> r != NULL)
        gsl_rng_free(mask -> r);
    if (mask -> blockMissing != NULL)
        free(mask -> blockMissing);
    if (mask -> correspondingSample != NULL)
        free(mask -> correspondingSample);
    if (mask -> prevCorrespondingSample != NULL)
        free(mask -> prevCorrespondingSample);
    if (kv_size(mask -> mask) != 0) {
        for (int i = 0; i < mask -> numRecords; i++)
            cb_destroy(kv_A(mask -> mask, i));
        kv_destroy(mask -> mask);
    }
    free(mask);
}