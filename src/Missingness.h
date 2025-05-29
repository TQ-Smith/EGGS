
// File: Missingness.h
// Date: 13 May 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Mask out genotypes in replicates.

// NOTE: I am iterating through the sample-by-record matrix a few more
//  times than I should. I keep it in an unoptimized form for now
//  to keep it straightforward.

#ifndef MISSINGNESS_H
#define MISSINGNESS_H

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#include "GenotypeParser.h"
#include "kissfft/kiss_fft.h"

// Holds the Fourier coefficients from a replicate.
typedef struct {
    int numSamples;
    int numRecords;
    kiss_fft_cpx** coeff;
} FourierCoefficients_t;

// NOTE: We can make Mask_t* an array of bit sets to increase efficiency.

// A Mask is a numRecords-by-numSamples matrix.
// An entry of 0 means not missing genotypes, and an entry of MISSING means
//  missing genotypes.
typedef struct {
    int numSamples;
    int numRecords;
    int** missing;
} Mask_t;

// Creates a genotype mask.
// Accepts:
//  int numSamples -> The number of samples.
//  int numRecords -> The number of records.
// Returns: Mask_t*, a numSamples-by-numRecords integer matrix.
Mask_t* init_mask(int numSamples, int numRecords);

// Calculate the Fourier coefficients for each sample's missing genotypes in a replicate.
FourierCoefficients_t* init_fourier_coefficients(Replicate_t* replicate, int numThreads);

// Create a numSamples-by-numRecords mask given a set of Fourier coefficients.
Mask_t* create_fourier_mask(FourierCoefficients_t* fourierCoeff, int numSamples, int numRecords, int numThreads);

// Create a random numSamples-by-numRecords Mask_t* from a beta distribution defined by mean and stder.
Mask_t* create_random_mask(int numSamples, int numRecords, double mean, double stder, int numThreads);

// Set loci between two adjacent missing loci to missing if base pair distance is <= fill.
void apply_fill(Replicate_t* replicate, Mask_t* mask, int fill);

// Free memory associated with mask.
void destroy_mask(Mask_t* mask);

// Free memory associated with fourierCoeff.
void destroy_fourier_coefficients(FourierCoefficients_t* fourierCoeff);

#endif