
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
#include <pthread.h>
#include "Sort.h"

Mask_t* init_mask(int numSamples, int numRecords) {
    Mask_t* mask = (Mask_t*) calloc(1, sizeof(Mask_t));
    mask -> numSamples = numSamples;
    mask -> numRecords = numRecords;
    mask -> missing = (int**) calloc(numSamples, sizeof(int*));
    for (int i = 0; i < numSamples; i++) {
        mask -> missing[i] = (int*) calloc(numRecords, sizeof(int));
    }
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

// Sort bins by power. Greatest to least. Quicksort.
void sort_bins(double* power, int* sortedBins, int numBins) {
    quick_sort_bins(power, sortedBins, 0, numBins - 1);
}

// Sort the first numBins from least to greatest.
void sort_indices(int* binIndices, int numBins) {
    quick_sort_indices(binIndices, 0, numBins - 1);
}

typedef struct {
    int startSampleIndex;
    int endSampleIndex;
    kiss_fft_scalar** inMask;
    FourierCoefficients_t* fourierCoeff;
} CoeffPackage_t;

void* fourier_coefficients(void* arg) {

    CoeffPackage_t* package = (CoeffPackage_t*) arg;

    kiss_fftr_cfg kiss_fft_state = kiss_fftr_alloc(package -> fourierCoeff -> numRecords, 0, 0, 0);
    // For each sample, forward transform the signal.
    for (int i = package -> startSampleIndex; i <= package -> endSampleIndex; i++) {
        package -> fourierCoeff -> coeff[i] = (kiss_fft_cpx*) calloc(package -> fourierCoeff -> numRecords / 2 + 1, sizeof(kiss_fft_cpx));
        kiss_fftr(kiss_fft_state, package -> inMask[i], package -> fourierCoeff -> coeff[i]);
    }

    free(kiss_fft_state);
    free(package);

    return NULL;
}

FourierCoefficients_t* init_fourier_coefficients(Replicate_t* replicate, int numThreads) {
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

    pthread_t* threads = NULL;
    int chunkSize = fourierCoeff -> numSamples / numThreads;
    int startSampleIndex = 0;

    // Assign each thread a group of samples.
    if (numThreads > 1) {
        threads = (pthread_t*) calloc(numThreads - 1, sizeof(pthread_t));
        for (int i = 0; i < numThreads - 1; i++) {
            CoeffPackage_t* package = calloc(1, sizeof(CoeffPackage_t));
            package -> inMask = inMask;
            package -> fourierCoeff = fourierCoeff;
            package -> startSampleIndex = startSampleIndex;
            package -> endSampleIndex = startSampleIndex + chunkSize - 1;
            pthread_create(&threads[i], NULL, fourier_coefficients, (void*) package);
            startSampleIndex += chunkSize;
        }
    }

    CoeffPackage_t* package = calloc(1, sizeof(CoeffPackage_t));
    package -> inMask = inMask;
    package -> fourierCoeff = fourierCoeff;
    package -> startSampleIndex = startSampleIndex;
    package -> endSampleIndex = fourierCoeff -> numSamples - 1;
    fourier_coefficients((void*) package);

    if (numThreads > 1) {
        for (int i = 0; i < numThreads - 1; i++)
            pthread_join(threads[i], NULL);
        free(threads);
    }

    for (int i = 0; i < fourierCoeff -> numSamples; i++)
        free(inMask[i]);
    free(inMask);

    return fourierCoeff;
}

typedef struct {
    int startSampleIndex;
    int endSampleIndex;
    int numSamples;
    int numRecords;
    int offset;
    FourierCoefficients_t* fourierCoeff;
    Mask_t* mask;
} FourierPackage_t;

#define POWER(coeff) (coeff.r * coeff.r + coeff.i * coeff.i)

void compressRep(kiss_fft_cpx* rep, double* power, int* powIndSorted, kiss_fft_cpx* freqs, gsl_rng* r, FourierPackage_t* package) {

    // Create our random replicate.
    double totalPower = 0;
    for (int i = 0; i < package -> fourierCoeff -> numRecords / 2 + 1; i++) {
        int randSample = (int) (package -> fourierCoeff -> numSamples * gsl_rng_uniform(r));
        rep[i].r = package -> fourierCoeff -> coeff[randSample][i].r;
        rep[i].i = package -> fourierCoeff -> coeff[randSample][i].i;
        power[i] = POWER(rep[i]) / package -> fourierCoeff -> numRecords;
        powIndSorted[i] = i;
        if (i == 0 || i == package -> fourierCoeff -> numRecords / 2)
            totalPower += power[i];
        else
            totalPower += 2 * power[i];
    }

    // Sort the bins and indices.
    sort_bins(power, powIndSorted, package -> fourierCoeff -> numRecords / 2 + 1);

    // Get the top bins in ascending order.
    sort_indices(powIndSorted, package -> numRecords / 2 + 1);

    // Calcuate compressed power.
    double compPower = 0;
    for (int i = 0; i < package -> numRecords / 2 + 1; i++)
        if (i == 0 || i == package -> numRecords / 2)
            compPower += POWER(rep[powIndSorted[i]]) / package -> numRecords;
        else
            compPower += 2 * POWER(rep[powIndSorted[i]]) / package -> numRecords;

    for (int i = 0; i < package -> numRecords / 2 + 1; i++) {
        freqs[i].r = rep[powIndSorted[i]].r;
        freqs[i].i = rep[powIndSorted[i]].i;
    }
}

void stretchRep(double* power, kiss_fft_cpx* freqs, gsl_rng* r, FourierPackage_t* package) {

    // Create our random replicate up to Nyquist.
    for (int i = 0; i < package -> fourierCoeff -> numRecords / 2 + 1; i++) {
        int randSample = (int) (package -> fourierCoeff -> numSamples * gsl_rng_uniform(r));
        freqs[i].r = package -> fourierCoeff -> coeff[randSample][i].r;
        freqs[i].i = package -> fourierCoeff -> coeff[randSample][i].i;
        power[i] = POWER(freqs[i]);
    }

    // Find the uniform tail of the power spectrum.

}

void* fourier_mask(void* arg) {

    FourierPackage_t* package = (FourierPackage_t*) arg;

    // Used to create the random replicate.
    kiss_fft_cpx* rep = (kiss_fft_cpx*) calloc(package -> fourierCoeff -> numRecords / 2 + 1, sizeof(kiss_fft_cpx));
    double* power = (double*) calloc(package -> fourierCoeff -> numRecords / 2 + 1, sizeof(double));
    int* powIndSorted = (int*) calloc(package -> fourierCoeff -> numRecords / 2 + 1, sizeof(int));

    // Used to stretch or compress replicate.
    kiss_fft_cpx* freqs = (kiss_fft_cpx*) calloc(package -> numRecords / 2 + 1, sizeof(kiss_fft_cpx));
    kiss_fft_scalar* inv = (kiss_fft_scalar*) calloc(package -> numRecords, sizeof(kiss_fft_scalar));
    kiss_fftr_cfg kiss_fft_state = kiss_fftr_alloc(package -> numRecords, 1, 0, 0);

    // Create our RNG.
    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));

    for (int i = package -> startSampleIndex; i <= package -> endSampleIndex; i++) {

        // If we are compressing or stretching the signal.
        if (package -> numRecords < package -> fourierCoeff -> numRecords)
            compressRep(rep, power, powIndSorted, freqs, r, package);
        else
            stretchRep(power, freqs, r, package);
        
        // Backward transfrom.
        kiss_fftri(kiss_fft_state, freqs, inv);
        // We do not need to normalize because the sign determines the state.
        for (int j = 0; j < package -> numRecords - (package -> offset); j++) {
            // If output signal is positive, then we set the mask element.
            if (inv[j] > 0)
                package -> mask -> missing[i][j] = MISSING;
            else 
                package -> mask -> missing[i][j] = 0;
        }
    }

    free(rep);
    free(power);
    free(powIndSorted);
    free(freqs);
    free(inv);
    free(kiss_fft_state);
    gsl_rng_free(r);
    free(package);

    return NULL;
}

Mask_t* create_fourier_mask(FourierCoefficients_t* fourierCoeff, int numSamples, int numRecords, int numThreads) {
    if (fourierCoeff == NULL || numSamples == 0 || numRecords == 0)
        return NULL;

    // Create an empty mask.
    Mask_t* mask = init_mask(numSamples, numRecords);

    // Same thing. Append empty record if odd number of records given.
    int offset = 1;
    if (numRecords % 2 == 0)
        offset = 0;
    numRecords += offset;

    pthread_t* threads = NULL;
    int chunkSize = numSamples / numThreads;
    int startSampleIndex = 0;

    // Assign each thread a group of samples.
    if (numThreads > 1) {
        threads = (pthread_t*) calloc(numThreads - 1, sizeof(pthread_t));
        for (int i = 0; i < numThreads - 1; i++) {
            FourierPackage_t* package = calloc(1, sizeof(FourierPackage_t));
            package -> numRecords = numRecords;
            package -> numSamples = numSamples;
            package -> offset = offset;
            package -> fourierCoeff = fourierCoeff;
            package -> mask = mask;
            package -> startSampleIndex = startSampleIndex;
            package -> endSampleIndex = startSampleIndex + chunkSize - 1;
            pthread_create(&threads[i], NULL, fourier_mask, (void*) package);
            startSampleIndex += chunkSize;
        }
    }

    FourierPackage_t* package = calloc(1, sizeof(FourierPackage_t));
    package -> numRecords = numRecords;
    package -> numSamples = numSamples;
    package -> offset = offset;
    package -> fourierCoeff = fourierCoeff;
    package -> mask = mask;
    package -> startSampleIndex = startSampleIndex;
    package -> endSampleIndex = numSamples - 1;
    fourier_mask((void*) package);

    if (numThreads > 1) {
        for (int i = 0; i < numThreads - 1; i++)
            pthread_join(threads[i], NULL);
        free(threads);
    }

    return mask;
}

typedef struct {
    double alpha;
    double beta;
    double* distribution;
    int sizeOfDistribution;
    int startRecordIndex;
    int endRecordIndex;
    Mask_t* mask;
} RandomPackage_t;

void* random_mask(void* arg) {

    RandomPackage_t* package = (RandomPackage_t*) arg;

    // Create our random number generator.
    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));

    // Create an array from 1 to numSamples.
    int* permu = (int*) calloc(package -> mask -> numSamples, sizeof(int));
    for (int i = 0; i < package -> mask -> numSamples; i++)
        permu[i] = i;

    // For each record in the partition.
    for (int i = package -> startRecordIndex; i <= package -> endRecordIndex; i++) {
        // Create random permutation of samples.
        shuffle_real_array(r, permu, package -> mask -> numSamples);
        // Draw a random proportion of missing samples from a the supplied distribution or beta distribution.
        int numMissing = 0; 
        if (package -> distribution == NULL)
            numMissing = (int) ((package -> mask -> numSamples) * gsl_ran_beta(r, package -> alpha, package -> beta));
        else
            numMissing = (int) ((package -> mask -> numSamples) * package -> distribution[(int) (package -> sizeOfDistribution * gsl_rng_uniform(r))]);
        // Set missing for samples at that record.
        for (int j = 0; j < numMissing; j++)
            package -> mask -> missing[permu[j]][i] = MISSING;
    }

    gsl_rng_free(r);
    free(permu);
    free(package);

    return NULL;
}

Mask_t* create_random_mask(int numSamples, int numRecords, double* distribution, int sizeOfDistribution, double mean, double stder, int numThreads) {

    // Create our mask.
    Mask_t* mask = init_mask(numSamples, numRecords);

    // Convert mean and stderr to alpha and beta that define a beta distribution.
    double alpha = -1;
    double beta = -1;
    if (distribution == NULL) {
        alpha = (mean * mean * (1 - mean)) / (stder * stder) - mean;
        beta = (alpha / mean) * (1 - mean);
    }

    // Unlike with the Fourier mask, we are partitioning w.r.t. records, not samples.

    pthread_t* threads = NULL;
    int chunkSize = numRecords / numThreads;
    int startRecordIndex = 0;

    // Assign each thread a partition of records.
    if (numThreads > 1) {
        threads = (pthread_t*) calloc(numThreads - 1, sizeof(pthread_t));
        for (int i = 0; i < numThreads - 1; i++) {
            RandomPackage_t* package = calloc(1, sizeof(RandomPackage_t));
            package -> alpha = alpha;
            package -> beta = beta;
            package -> distribution = distribution;
            package -> sizeOfDistribution = sizeOfDistribution;
            package -> mask = mask;
            package -> startRecordIndex = startRecordIndex;
            package -> endRecordIndex = startRecordIndex + chunkSize - 1;
            pthread_create(&threads[i], NULL, random_mask, (void*) package);
            startRecordIndex += chunkSize;
        }
    }

    RandomPackage_t* package = calloc(1, sizeof(RandomPackage_t));
    package -> alpha = alpha;
    package -> beta = beta;
    package -> distribution = distribution;
    package -> sizeOfDistribution = sizeOfDistribution;
    package -> mask = mask;
    package -> startRecordIndex = startRecordIndex;
    package -> endRecordIndex = numRecords - 1;
    random_mask((void*) package);

    if (numThreads > 1) {
        for (int i = 0; i < numThreads - 1; i++)
            pthread_join(threads[i], NULL);
        free(threads);
    }

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
            if (mask -> missing[i][j] == MISSING) {
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
        for (int i = 0; i < mask -> numSamples; i++)
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