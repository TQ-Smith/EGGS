
// File: Missingness.c
// Date: 13 May 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Mask out genotypes in replicates.

#include "Missingness.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include <time.h>
#include "kissfft/kiss_fftr.h"

Mask_t* init_mask(int numSamples, int numRecords) {
    Mask_t* mask = (Mask_t*) calloc(1, sizeof(Mask_t));
    mask -> numSamples = numSamples;
    mask -> numRecords = numRecords;
    mask -> missing = (int**) calloc(numRecords, sizeof(int*));
    for (int i = 0; i < numRecords; i++) {
        mask -> missing[i] = (int*) calloc(numSamples, sizeof(int));
    }
    return mask;
}

// Fisher-Yates shuffle algorithm for an integer array.
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
    if (replicate == NULL)
        return NULL;
    if (replicate -> numRecords == 0 || replicate -> numSamples == 0)
        return NULL;

    // If the replicate has an odd number of records, then we append a record without missing genotypes.
    //  We do this because the FFT for a real sequence requires an even number of elements.
    int offset = 1;
    if ((replicate -> numRecords) % 2 == 0)
        offset = 0;

    // For each sample, we create an array spanning the number of records.
    kiss_fft_scalar** inMask = (kiss_fft_scalar**) calloc(replicate -> numSamples, sizeof(kiss_fft_scalar*));
    for (int i = 0; i < replicate -> numSamples; i++) {
        inMask[i] = (kiss_fft_scalar*) calloc(replicate -> numRecords + offset, sizeof(kiss_fft_scalar));
        int j = 0;
        for (Record_t* head = replicate -> headRecord; head != NULL; head = head -> nextRecord) {
            // A missing entry is given a 1, and a non-missing entry is given a -1.
            if (head -> genotypes[i].left == MISSING && head -> genotypes[i].right == MISSING) 
                inMask[i][j] = 1;
            else
                inMask[i][j] = -1;
            j++;
        }
        // If a dummy record was appended, then set all genotypes to non-missing.
        if (offset == 1)
            inMask[i][replicate -> numRecords] = -1;
    }

    // Create our set of Fourier coefficients.
    FourierCoefficients_t* fourierCoeff = (FourierCoefficients_t*) calloc(1, sizeof(FourierCoefficients_t));
    fourierCoeff -> numSamples = replicate -> numSamples;
    fourierCoeff -> numRecords = replicate -> numRecords + offset;
    fourierCoeff -> coeff = (kiss_fft_cpx**) calloc(fourierCoeff -> numSamples, sizeof(kiss_fft_cpx*));
    kiss_fftr_cfg kiss_fft_state = kiss_fftr_alloc(fourierCoeff -> numRecords, 0, 0, 0);

    // For each sample, forward transform the signal.
    for (int i = 0; i < fourierCoeff -> numSamples; i++) {
        fourierCoeff -> coeff[i] = (kiss_fft_cpx*) calloc(fourierCoeff -> numRecords / 2 + 1, sizeof(kiss_fft_cpx));
        kiss_fftr(kiss_fft_state, inMask[i], fourierCoeff -> coeff[i]);
    }

    free(kiss_fft_state);
    for (int i = 0; i < fourierCoeff -> numSamples; i++)
        free(inMask[i]);
    free(inMask);
    kiss_fft_cleanup();

    return fourierCoeff;
}

Mask_t* create_fourier_mask(FourierCoefficients_t* fourierCoeff, int numSamples, int numRecords) {
    if (fourierCoeff == NULL || numSamples == 0 || numRecords == 0)
        return NULL;

    // Create an empty mask.
    Mask_t* mask = init_mask(numSamples, numRecords);

    // Same thing. Append empty record if odd number of records given.
    int offset = 1;
    if (numRecords % 2 == 0)
        offset = 0;
    numRecords += offset;

    kiss_fft_cpx* freqs = (kiss_fft_cpx*) calloc(numRecords / 2 + 1, sizeof(kiss_fft_cpx));
    kiss_fft_scalar* inv = (kiss_fft_scalar*) calloc(numRecords, sizeof(kiss_fft_scalar));
    kiss_fftr_cfg kiss_fft_state = kiss_fftr_alloc(numRecords, 1, 0, 0);

    // Create our RNG.
    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));

    // We preserve amplitude on backward transform. This is not needed since the sign of the
    //  element determines the state.
    double normalizationFactor = 1.0 / fourierCoeff -> numRecords;
    if (numRecords < fourierCoeff -> numRecords)
        normalizationFactor = 1.0 / numRecords;

    for (int i = 0; i < numSamples; i++) {
        for (int j = 0; j < numRecords / 2 + 1; j++) {
            // Truncate spectrum if output signal is shorter then input signal.
            if (j <= fourierCoeff -> numRecords / 2) {
                // Randomly choose sample's Fourier coefficient.
                int randSample = (int) (fourierCoeff -> numSamples * gsl_rng_uniform(r));
                freqs[j].r = fourierCoeff -> coeff[randSample][j].r;
                freqs[j].i = fourierCoeff -> coeff[randSample][j].i;
            // Zero pad if output signal is longer than input signal.
            } else {
                freqs[j].r = 0;
                freqs[j].i = 0;
            }
        }
        // printf("Freqs:\n");
        // for (int j = 0; j < numRecords / 2 + 1; j++) {
           // printf("%.3f+%.3fi\n", freqs[j].r, freqs[j].i);
        // }

        // Backward transfrom.
        kiss_fftri(kiss_fft_state, freqs, inv);
        // printf("Missing:\n");
        for (int j = 0; j < numRecords - offset; j++) {
            inv[j] *= normalizationFactor;
            // If output signal is positive, then we set the mask element.
            if (inv[j] > 0)
                mask -> missing[j][i] = MISSING;
            else 
                mask -> missing[j][i] = 0;
            // printf("%.3f\t%d\n", inv[j], mask -> missing[j][i]);
        }
    }

    free(freqs);
    free(inv);
    free(kiss_fft_state);
    gsl_rng_free(r);
    kiss_fft_cleanup();
    return mask;
}

Mask_t* create_random_mask(int numSamples, int numRecords, double mean, double stder) {
    // Create our random number generator.
    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));

    // Create our mask.
    Mask_t* mask = init_mask(numSamples, numRecords);

    // Create an array from 1 to numRecords.
    int* permu = (int*) calloc(numRecords, sizeof(int));
    for (int i = 0; i < numRecords; i++)
        permu[i] = i + 1;

    // COnvert mean and stderr to alpha and beta that define a beta distribution.
    double alpha = (mean * mean * (1 - mean)) / (stder * stder) - mean;
    double beta = (alpha / mean) * (1 - mean);

    // For each sample.
    for (int i = 0; i < numSamples; i++) {
        // Create random permutation of records.
        shuffle_real_array(r, permu, numRecords);
        // Draw a random proportion of missing records from a beta distribution.
        int numMissing = (int) (numRecords * gsl_ran_beta(r, alpha, beta));
        // Set records to missing for that sample.
        for (int j = 0; j < numMissing; j++)
            mask -> missing[permu[j] - 1][i] = MISSING;
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

void destroy_fourier_coefficients(FourierCoefficients_t* fourierCoeff) {
    if (fourierCoeff == NULL)
        return;
    if (fourierCoeff -> coeff != NULL) {
        for (int i = 0; i < fourierCoeff -> numSamples; i++)
            free(fourierCoeff -> coeff[i]);
        free(fourierCoeff -> coeff);
    }
    free(fourierCoeff);
}