
#include "GenotypeParser.h"

InputStream_t* init_input_stream(FILE* source) {
    if (source == NULL)
        return NULL;
    InputStream_t* inputStream = calloc(1, sizeof(InputStream_t));
    int fIn = fileno(source);
    inputStream -> file = gzdopen(fIn, "r");
    inputStream -> fpIn = ks_init(inputStream -> file);
    inputStream -> buffer = calloc(1, sizeof(kstring_t));
    return inputStream;
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
    free(inputStream);
}

Replicate_t* init_vcf_replicate(InputStream_t* inputStream) {
    if (inputStream == NULL)
        return NULL;

    Replicate_t* replicate = calloc(1, sizeof(Replicate_t));

    int dret;
    do {
        ks_getuntil(inputStream -> fpIn, '\n', inputStream -> buffer, &dret);
    } while (strncmp(inputStream -> buffer -> s, "#C", 2) != 0);

    int numSamples = 0;
    for (int i = 0; i < inputStream -> buffer -> l; i++)
        if (inputStream -> buffer -> s[i] == '\t')
            numSamples++;
    replicate -> numSamples = numSamples - 8;


    int numTabs = 0, prevIndex = 0;
    replicate -> sampleNames = calloc(replicate -> numSamples, sizeof(char*));
    for (int i = 0; i <= inputStream -> buffer -> l; i++) {
        if (i == inputStream -> buffer -> l || inputStream -> buffer -> s[i] == '\t') {
            if (numTabs > 8) {
                replicate -> sampleNames[numTabs - 9] = strndup(inputStream -> buffer -> s + prevIndex + 1, i - prevIndex - 1);
            }
            prevIndex = i;
            numTabs++;
        }
    }

    return replicate;
}

bool get_next_vcf_record(Record_t* record, InputStream_t* inputStream) {

    int numTabs = 0, prevIndex = 0, dret = 0, numAlleles = 2;

    ks_getuntil(inputStream -> fpIn, '\n', inputStream -> buffer, &dret);

    if (ks_eof(inputStream -> fpIn) || inputStream -> buffer -> l == 0) 
        return false;

    for (int i = 0; i <= inputStream -> buffer -> l; i++) {
        if (i == inputStream -> buffer -> l || inputStream -> buffer -> s[i] == '\t') {
            if (numTabs == 0) {
                if (record -> chrom != NULL)
                    free(record -> chrom);
                record -> chrom = strndup(inputStream -> buffer -> s, i);
            } else if (numTabs == 1) {
                record -> position = (int) strtol(inputStream -> buffer -> s + prevIndex + 1, (char**) NULL, 10);
            } else if (numTabs == 3) {
                if (record -> ref != NULL)
                    free(record -> ref);
                record -> ref = strndup(inputStream -> buffer -> s + prevIndex + 1, i - prevIndex - 1);
            } else if (numTabs == 4) {
                for (int j = prevIndex + 1; inputStream -> buffer -> s[j] != '\t'; j++)
                    if (inputStream -> buffer -> s[j] == ',')
                        numAlleles++;
                if (record -> alts != NULL)
                    free(record -> alts);
                record -> alts = strndup(inputStream -> buffer -> s + prevIndex + 1, i - prevIndex - 1);
            } else if (numTabs > 8) {
                record -> genotypes[numTabs - 9].left = -1;
                record -> genotypes[numTabs - 9].right = -1;
                record -> genotypes[numTabs - 9].isPhased = false;
                char* start = inputStream -> buffer -> s + prevIndex + 1;
                char* next = start + 1;
                if (start[0] != '.')
                    record -> genotypes[numTabs - 9].left = (int) strtol(start, &next, 10);
                if (next[0] == '|')
                    record -> genotypes[numTabs - 9].isPhased = true;
                if ((next[0] == '|' || next[0] == '/') && next[1] != '.')
                    record -> genotypes[numTabs - 9].right = (int) strtol(next + 1, (char**) NULL, 10);
            }
            prevIndex = i;
            numTabs++;
        }
    }
    record -> numAlleles = numAlleles;

    return true;

}

void parse_vcf(Replicate_t* replicate, InputStream_t* inputStream) {

    while (true) {

        Record_t* record = calloc(1, sizeof(Record_t));
        record -> genotypes = calloc(replicate -> numSamples, sizeof(Genotype_t));

        if (get_next_vcf_record(record, inputStream)) {
            if (replicate -> numRecords == 0) {
                replicate -> headRecord = record;
                replicate -> tailRecord = record;
            } else {
                replicate -> tailRecord -> nextRecord = record;
                replicate -> tailRecord = record;
            }
            replicate -> numRecords++;
        } else {
            free(record -> genotypes);
            free(record);
            break;
        }
        
    }

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
        replicate -> headRecord = replicate -> headRecord -> nextRecord;
        destroy_record(temp);
    }
    if (replicate -> sampleNames != NULL) {
        for (int i = 0; i < replicate -> numSamples; i++) {
            free(replicate -> sampleNames[i]);       
        }
        free(replicate -> sampleNames);
    }
    free(replicate);
}