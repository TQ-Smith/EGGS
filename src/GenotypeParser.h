
// File: GenotypeParser.h
// Date: 13 May 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Parse VCF files and ms-style replicates.

#ifndef GENOTYPE_PARSER_H
#define GENOTYPE_PARSER_H

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#include <stdbool.h>
#include "kseq.h"
#include <zlib.h>
#include <stdio.h>

// We use -1 to denote a missing genotype.
#define MISSING -1

// A stream to read GZIP files.
#define BUFFER_SIZE 4096
KSTREAM_INIT(gzFile, gzread, BUFFER_SIZE)

// Represents the genotype of a sample at a locus.
typedef struct {
    // If the genotype is phased.
    bool isPhased;
    // The left and right genotypes.
    int left;
    int right;
} Genotype_t;

// Represents a record/locus.
typedef struct Record {
    // The chromosome.
    char* chrom;
    // The position on the chromosome.
    int position;
    // The reference allele.
    char* ref;
    // The alternative alleles. Seperated by ',' is multiple.
    char* alts;
    // The total number of alleles at the locus.
    int numAlleles;
    // Redunant but useful. The number of samples in the genotypes array.
    int numSamples;
    // The index of the record in the list.
    int recordIndex;
    // The genotypes for each sample.
    Genotype_t* genotypes;
    // A pointer to the next locus. Record_t is a node in a linked list.
    struct Record* nextRecord;
} Record_t;


// A replicate is a linked list of records.
typedef struct {
    int numSamples;
    // Number of records in the list.
    int numRecords;
    // The sample names.
    char** sampleNames;
    // Head and tail pointers.
    Record_t* headRecord;
    Record_t* tailRecord;
} Replicate_t;

// Packages all the necessary information to read from a GZIP file.
typedef struct {
    // The file we are reading.
    gzFile file;
    // The stream that wraps the file.
    kstream_t* fpIn;
    // The buffer to read in lines of the file.
    kstring_t* buffer;
} InputStream_t;

// From a file name, open up a GZIP stream.
// Accepts:
//  char* source -> The file name open. If null, read from stdin.
// Returns: InputStream_t*, the input stream we are reading from. NULL, if source does not exist. 
InputStream_t* init_input_stream(char* source);

// Closes file stream and deallocates memory.
void destroy_input_stream(InputStream_t* inputStream);

// From an input stream reading a VCF file, create an empty replicate.
// Accepts:
//  InputStream_t* inputStream -> The input stream to read from.
// Returns: Replicate_t*, an empty replicate with sampleNames set from VCF header info.
Replicate_t* init_vcf_replicate(InputStream_t* inputStream);

// Read in next VCF record from input stream.
// Accepts:
//  Record_t* record -> An allocated record. Set fields.
//  InputStream_t* inputStream -> The VCF file we are reading.
// Returns: bool, true is a record was read. false, if eof.
bool get_next_vcf_record(Record_t* record, InputStream_t* inputStream);

// Reads in whole VCF file and adds each record to replicate.
// Accepts:
//  Replicate_t* replicate -> An allocated, but empty, replicate to add records to.
//  InputStream_t* inputStream -> The VCF file we are reading from.
// Returns: void.
void parse_vcf(Replicate_t* replicate, InputStream_t* inputStream);

// Reads in a ms replicate from input stream and returns replicate.
// Accepts:
//  InputStream_t* inputStream -> The ms-style file we are reading from.
//  int length -> The length of the resulting segment in base pairs.
// Returns: Replicate_t*, the ms replicate. NULL, if no more ms replicates in stream. 
Replicate_t* parse_ms(InputStream_t* inputStream, int length);

// Frees memory associated with a record.
void destroy_record(Record_t* record);

// Frees memory associated with a replicate.
void destroy_replicate(Replicate_t* replicate);

#endif