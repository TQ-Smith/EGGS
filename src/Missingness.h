
#ifndef MISSINGNESS_H
#define MISSINGNESS_H

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#include "GenotypeParser.h"
#include "kissfft/kiss_fft.h"

typedef struct {
    int numSamples;
    int numRecords;
    kiss_fft_cpx** coeff;
} FourierCoefficients_t;

typedef struct {
    int numSamples;
    int numRecords;
    int** missing;
} Mask_t;

Mask_t* init_mask(int numSamples, int numRecords);

FourierCoefficients_t* init_fourier_coefficients(Replicate_t* replicate);

Mask_t* create_fourier_mask(FourierCoefficients_t* fourierCoeff, int numSamples, int numRecords);

Mask_t* create_random_mask(int numSamples, int numRecords, double mean, double stder);

void apply_fill(Replicate_t* replicate, Mask_t* mask, int fill);

void destroy_mask(Mask_t* mask);

void destroy_fourier_coefficients(FourierCoefficients_t* fourierCoeff);

#endif