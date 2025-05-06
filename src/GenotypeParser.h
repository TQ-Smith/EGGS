#ifndef GENOTYPE_PARSER_H
#define GENOTYPE_PARSER_H

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#include <stdbool.h>
#include "kseq.h"
#include <zlib.h>
#include <stdio.h>

#define MISSING -1

#define BUFFER_SIZE 4096
KSTREAM_INIT(gzFile, gzread, BUFFER_SIZE)

typedef struct {
    bool isPhased;
    int left;
    int right;
} Genotype_t;

typedef struct Record {
    char* chrom;
    int position;
    char* ref;
    char* alts;
    int numAlleles;
    int numSamples;
    Genotype_t* genotypes;
    struct Record* nextRecord;
} Record_t;

typedef struct {
    int numSamples;
    int numRecords;
    char** sampleNames;
    Record_t* headRecord;
    Record_t* tailRecord;
} Replicate_t;

typedef struct {
    gzFile file;
    kstream_t* fpIn;
    kstring_t* buffer;
} InputStream_t;

InputStream_t* init_input_stream(FILE* source);

void destroy_input_stream(InputStream_t* inputStream);

Replicate_t* init_vcf_replicate(InputStream_t* inputStream);

bool get_next_vcf_record(Record_t* record, InputStream_t* inputStream);

void parse_vcf(Replicate_t* replicate, InputStream_t* inputStream);

Replicate_t* parse_ms(InputStream_t* inputStream, int length);

void destroy_record(Record_t* record);

void destroy_replicate(Replicate_t* replicate);

#endif