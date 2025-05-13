
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

    int offset = 1;
    if ((replicate -> numRecords) % 2 == 0)
        offset = 0;

    kiss_fft_scalar** inMask = (kiss_fft_scalar**) calloc(replicate -> numSamples, sizeof(kiss_fft_scalar*));
    for (int i = 0; i < replicate -> numSamples; i++) {
        inMask[i] = (kiss_fft_scalar*) calloc(replicate -> numRecords + offset, sizeof(kiss_fft_scalar));
        int j = 0;
        for (Record_t* head = replicate -> headRecord; head != NULL; head = head -> nextRecord) {
            if (head -> genotypes[i].left == MISSING && head -> genotypes[i].right == MISSING) 
                inMask[i][j] = 1;
            else
                inMask[i][j] = -1;
            j++;
        }
        if (offset == 1)
            inMask[i][replicate -> numRecords] = -1;
    }

    FourierCoefficients_t* fourierCoeff = (FourierCoefficients_t*) calloc(1, sizeof(FourierCoefficients_t));
    fourierCoeff -> numSamples = replicate -> numSamples;
    fourierCoeff -> numRecords = replicate -> numRecords + offset;
    fourierCoeff -> coeff = (kiss_fft_cpx**) calloc(fourierCoeff -> numSamples, sizeof(kiss_fft_cpx*));
    kiss_fftr_cfg kiss_fft_state = kiss_fftr_alloc(fourierCoeff -> numRecords, 0, 0, 0);

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

    Mask_t* mask = init_mask(numSamples, numRecords);

    int offset = 1;
    if (numRecords % 2 == 0)
        offset = 0;
    numRecords += offset;

    kiss_fft_cpx* freqs = (kiss_fft_cpx*) calloc(numRecords / 2 + 1, sizeof(kiss_fft_cpx));
    kiss_fft_scalar* inv = (kiss_fft_scalar*) calloc(numRecords, sizeof(kiss_fft_scalar));
    kiss_fftr_cfg kiss_fft_state = kiss_fftr_alloc(numRecords, 1, 0, 0);

    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));

    double normalizationFactor = 1.0 / fourierCoeff -> numRecords;
    if (numRecords < fourierCoeff -> numRecords)
        normalizationFactor = 1.0 / numRecords;

    for (int i = 0; i < numSamples; i++) {
        for (int j = 0; j < numRecords / 2 + 1; j++) {
            if (j <= fourierCoeff -> numRecords / 2) {
                int randSample = (int) (fourierCoeff -> numSamples * gsl_rng_uniform(r));
                freqs[j].r = fourierCoeff -> coeff[randSample][j].r;
                freqs[j].i = fourierCoeff -> coeff[randSample][j].i;
            } else {
                freqs[j].r = 0;
                freqs[j].i = 0;
            }
        }
        // printf("Freqs:\n");
        // for (int j = 0; j < numRecords / 2 + 1; j++) {
           // printf("%.3f+%.3fi\n", freqs[j].r, freqs[j].i);
        // }
        kiss_fftri(kiss_fft_state, freqs, inv);
        // printf("Missing:\n");
        for (int j = 0; j < numRecords - offset; j++) {
            inv[j] *= normalizationFactor;
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
    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));

    Mask_t* mask = init_mask(numSamples, numRecords);

    int* permu = (int*) calloc(numRecords, sizeof(int));
    for (int i = 0; i < numRecords; i++)
        permu[i] = i + 1;

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
    Record_t* temp = NULL;
    int prevMisGenoPos = -1;
    for (int i = 0; i < mask -> numSamples; i++) { 
        temp = replicate -> headRecord;
        prevMisRecord = NULL;
        for (int j = 0; j < mask -> numRecords; j++) {
            if (mask -> missing[j][i] == MISSING) {
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