
// File: Missingness.h
// Date: 13 May 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Mask out genotypes in replicates.

#ifndef MISSINGNESS_H
#define MISSINGNESS_H

#include "compact_bitset.h"
#include "GenotypeParser.h"
#include "kvec.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

// Our mask structure.
typedef struct {
    // For each record, set bit if sample is missing.
    kvec_t(CompactBitset*) mask;
    // Number of samples/records in the VCF.
    int numSamples;
    int numRecords;
    gsl_rng* r;  
    // Holds by sample prob of missing.
    double* blockMissing;
    int* correspondingSample;
    int* prevCorrespondingSample;
} MissingMask_t;

// Fisher-Yates shuffle algorithm for an integer array.
void shuffle_real_array(gsl_rng* r, int* array, int n);

// Create our mask.
// Accepts:
//  Replicate_t* replicate -> A replicate. Only header is read.
//  InputStream_t* inputStream -> The input stream to the VCF.
// Returns: MissingMask_t*, sets bit for sample if missing at each site.
MissingMask_t* init_missing_mask(Replicate_t* replicate, InputStream_t* inputStream);

// Get the mask for the current site.
// Accepts:
//  MissingMask_t* mask -> Our mask to replicate.
//  CompactBitset* cb -> Sets bit for samples to mark as missing in current record.
//  int numRecords -> numRecords for the target.
//  int numSaamples -> numSamples for the target.
//  int site -> The current record.
void get_mask_for_next_site(MissingMask_t* mask, CompactBitset* cb, int numRecords, int numSamples, int site);

// Free the mask structure.
void destroy_missing_mask(MissingMask_t* mask);

#endif