
#include "Missingness.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include <time.h>

Mask_t* init_mask(int numSamples, int numRecords) {
    Mask_t* mask = calloc(1, sizeof(Mask_t));
    mask -> numSamples = numSamples;
    mask -> numRecords = numRecords;
    mask -> missing = calloc(numRecords, sizeof(int*));
    for (int i = 0; i < numRecords; i++) {
        mask -> missing[i] = calloc(numSamples, sizeof(int));
    }
    return mask;
}

void shuffle_real_array(gsl_rng* r, int* array, int n) {
    int temp = 0, j = 0;
    for (int i = 0; i < n - 1; i++) {
        j = (int) (i + gsl_rng_uniform(r) / (RAND_MAX / (n - i) + 1));
        temp = array[j];
        array[j] = array[i];
        array[i] = temp;
    }
}

FourierCoefficients_t* init_fourier_coefficients(Replicate_t* replicate) {
    return NULL;
}

Mask_t* create_fourier_mask(FourierCoefficients_t* fourierCoeff, int numSamples, int numRecords) {
    return NULL;
}

Mask_t* create_random_mask(int numSamples, int numRecords, double mean, double stder) {
    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));

    Mask_t* mask = init_mask(numSamples, numRecords);

    int* permu = calloc(numRecords, sizeof(int));
    for (int i = 1; i <= numRecords; i++)
        permu[i] = i;

    double alpha = (mean * mean * (1 - mean)) / (stder * stder) - mean;
    double beta = (alpha / mean) * (1 - mean);

    for (int i = 0; i < numSamples; i++) {
        shuffle_real_array(r, permu, numRecords);
        int numMissing = (int) (numRecords * gsl_ran_beta(r, alpha, beta));
        for (int j = 0; j < numMissing; j++)
            mask -> missing[permu[j]][i] = MISSING;
    }

    gsl_rng_free(r);
    free(permu);

    return mask;
}

void apply_fill(Replicate_t* replicate, Mask_t* mask, int fill) {
    Record_t* prevMisRecord = NULL;
    Record_t* temp = replicate -> headRecord;
    int prevMisGenoPos = -1;
    for (int i = 0; i < mask -> numSamples; i++) { 
        for (int j = 0; j < mask -> numRecords; j++) {
            if (prevMisGenoPos == -1 || mask -> missing[j][i] == MISSING) {
                if (prevMisGenoPos != -1 && (temp -> position - prevMisGenoPos) <= fill) {
                    for (Record_t* k = prevMisRecord -> nextRecord; k != temp; k = k -> nextRecord) {
                        k -> genotypes[j].left = MISSING;
                        k -> genotypes[j].right = MISSING;
                    }
                } else {
                    prevMisGenoPos = temp -> position;
                    prevMisRecord = temp;
                }
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

void destroy_fourier_coefficients(int numSamples) {

}