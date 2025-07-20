
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
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

// The proportion of missing samples at each site.
typedef struct {
    double* proportions;
    int numRecords;
    gsl_rng* r;
    int* permu;
    int sizeOfPermu;
} MissingDistribution_t;

// Calculate the proportion of missing samples at each site.
// Accepts:
//  Recplicate_t* replicate -> A replicate. Only header is read.
//  InputStream_t* inputStream -> The input stream to the VCF.
// Returns: MissingDistribution_t*, the proportion of samples missing at each site.
MissingDistribution_t* init_missing_distribution(Replicate_t* replicate, InputStream_t* inputStream);

// Create a mask for a locus.
// Accepts:
//  MissingDistribution_t* dis -> The distribution used for random/EGGS mask algorithm.
//  int numSamples -> The number of samples at the locus.
//  int curRecord -> Used for EGGS mask algorithm. The current record the mask will be applied to.
//  int numRecords -> The total number of records.
//  double mean -> If beta distribution, then the mean.
//  double stder -> If beta distribution, then the stder.
//  int* mask -> Array of size numSamples. Element set to 1 if the sample is missing.
// Returns: void.
void create_random_mask(MissingDistribution_t* dis, int numSamples, int curRecord, int numRecords, double mean, double stder, int* mask);

// Free the distribution structure.
void destroy_missing_distribution(MissingDistribution_t* dis);

#endif