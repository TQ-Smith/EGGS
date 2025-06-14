
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

#include "GenotypeParser.h"

// NOTE: We can make Mask_t* an array of bit sets to increase efficiency.

// The proportion of missing samples at each site.
typedef struct {
    double* proportions;
    int numRecords;
} MissingDistribution_t;

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
// Returns: Mask_t*, a numRecords-by-numSamples integer matrix.
Mask_t* init_mask(int numSamples, int numRecords);

// Calculate the proportion of missing samples at each site.
// Accepts:
//  Recplicate_t* replicate -> A replicate read into memory.
// Returns: MissingDistribution_t*, the proportion of samples missing at each site.
MissingDistribution_t* init_missing_distribution(Replicate_t* replicate);

// Create a mask by replicating the structure of missingness.
// Accepts:
//  MissingDistribution_t* dis -> The proportion of missingness across samples.
//  int numSamples -> The number of samples in the mask.
//  int numRecords -> The number of records in the mask.
// Returns: Mask_t*, the created mask.
Mask_t* create_missing_mask(MissingDistribution_t* dis, int numSamples, int numRecords);

// Randomly create a mask from a predefined distribution or from the beta distribution defined by mean and stder.
Mask_t* create_random_mask(MissingDistribution_t* dis, int numSamples, int numRecords, double mean, double stder);

// Free memory associated with mask.
void destroy_mask(Mask_t* mask);

// Free the distribution structure.
void destroy_missing_distribution(MissingDistribution_t* dis);

#endif