
#include "GenotypeParser.h"

InputStream_t* init_input_stream(FILE* source) {
    InputStream_t* inputStream = calloc(1, sizeof(InputStream_t));
    if (source == NULL) {
        int fIn = fileno(stdin);
        inputStream -> file = gzdopen(fIn, "r");
        inputStream -> fpIn = ks_init(inputStream -> file);
        inputStream -> buffer = calloc(1, sizeof(kstring_t));
        return inputStream;
    } else { return NULL; }
}

void destroy_input_stream(InputStream_t* inputStream) {
    if (inputStream == NULL)
        return;
    if (inputStream -> fpIn != NULL)
        ks_destroy(inputStream -> fpIn);
    if (inputStream -> buffer != NULL) {
        if (inputStream -> buffer -> s != NULL)
            free(inputStream -> buffer -> s);
        free(inputStream -> buffer);
    }
    gzclose(inputStream -> file);
}

void get_next_vcf_record(Record_t* record, InputStream_t* inputStream) {

}

Replicate_t* parse_vcf(InputStream_t* inputStream) {
    if (inputStream == NULL)
        return NULL;
    return NULL;
}

Replicate_t* parse_ms(InputStream_t* inputStream) {
    if (inputStream == NULL)
        return NULL;
    return NULL;
}

void destroy_record(Record_t* record) {
    if (record == NULL)
        return;
    if (record -> genotypes != NULL)
        free(record -> genotypes);
    if (record -> chrom != NULL)
        free(record -> chrom);
    if (record -> ref != NULL)
        free(record -> ref);
    if (record -> alts != NULL)
        free(record -> alts);
    free(record);
}

void destroy_replicate(Replicate_t* replicate) {
    if (replicate == NULL)
        return;
    Record_t* temp = NULL;
    for (int i = 0; i < replicate -> numRecords; i++) {
        temp = replicate -> headRecord;
        replicate -> headRecord = replicate -> headRecord -> next_record;
        destroy_record(temp);
    }
    if (replicate -> sampleNames != NULL) {
        for (int i = 0; i < replicate -> numSamples; i++) {
            free(replicate -> sampleNames[i]);       
        }
        free(replicate -> sampleNames);
    }
}