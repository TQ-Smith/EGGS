
#ifndef MISSINGNESS_H
#define MISSINGNESS_H

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#include "GenotypeParser.h"
#include "fftpack/fft.h"

typedef struct {
    int numSamples;
    int numLoci;
    FFTTransformer** transformers;
} FourierCoefficients_t;

typedef struct {
    int numSamples;
    int numLoci;
    int** missing;
} Mask_t;

FourierCoefficients_t* init_fourier_coefficients(Replicate_t* replicate, int numThreads);

Mask_t* create_fourier_mask(FourierCoefficients_t* fourierCoeff, int numSamples, int numLoci, int numThreads);

Mask_t* create_random_mask(int numSamples, int numLoci, int numThreads);

void apply_fill(Replicate_t* replicate, Mask_t* mask, int fill);

void destroy_mask(Mask_t* mask);

void destroy_fourier_coefficients(int numSamples);

#endif