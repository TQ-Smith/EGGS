
// File: GenotypeParser.h
// Date: 13 May 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Parse VCF files and ms-style replicates.

#include "GenotypeParser.h"
#include "kvec.h"

InputStream_t* init_input_stream(char* source) {
    InputStream_t* inputStream = (InputStream_t*) calloc(1, sizeof(InputStream_t));
    // If source is NULL, read from stdin.
    if (source == NULL) {
        int fIn = fileno(stdin);
        inputStream -> file = gzdopen(fIn, "r");
    // Otherwise, we are reading from a file.
    } else 
        inputStream -> file = gzopen(source, "r");
    inputStream -> fpIn = ks_init(inputStream -> file);
    inputStream -> buffer = (kstring_t*) calloc(1, sizeof(kstring_t));
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

    // Allocate our replicate.
    Replicate_t* replicate = (Replicate_t*) calloc(1, sizeof(Replicate_t));

    // Eat meta info lines until line starting with #CHROM.
    int dret;
    do {
        ks_getuntil(inputStream -> fpIn, '\n', inputStream -> buffer, &dret);
    } while (strncmp(inputStream -> buffer -> s, "#C", 2) != 0);

    // Count the number of samples in the VCF file.
    int numSamples = 0;
    for (int i = 0; i < inputStream -> buffer -> l; i++)
        if (inputStream -> buffer -> s[i] == '\t')
            numSamples++;
    replicate -> numSamples = numSamples - 8;

    // Copy the sample names from the header line.
    replicate -> sampleNames = calloc(replicate -> numSamples, sizeof(char*));
    // Read in the sample names.
    char* header = strdup(inputStream -> buffer -> s);
    char* tok = NULL;
    tok = strtok(header, "\t");
    for (int i = 0; i < 9; i++)
        tok = strtok(NULL, "\t");
    for (int i = 0; i < replicate -> numSamples; i++) {
        replicate -> sampleNames[i] = strdup(tok);
        tok = strtok(NULL, "\t");
    }

    return replicate;
}

bool get_next_vcf_record(Record_t* record, InputStream_t* inputStream) {

    int dret = 0, numAlleles = 2;

    // Get the next line.
    ks_getuntil(inputStream -> fpIn, '\n', inputStream -> buffer, &dret);

    // If end of file, return false.
    if (ks_eof(inputStream -> fpIn) || inputStream -> buffer -> l == 0) 
        return false;

    // We parse the record. This is messier than splitting on tab.
    char* line = strdup(inputStream -> buffer -> s);
    char* tok = strtok(line, "\t");
    if (record -> chrom != NULL)
        free(record -> chrom);
    record -> chrom = strdup(tok);
    for (int i = 1; i < 9 + record -> numSamples; i++) {
        tok = strtok(NULL, "\t");
        if (i == 1)
            record -> position = (int) strtol(tok, (char**) NULL, 10);
        else if (i == 3) {
            if (record -> ref != NULL)
                free(record -> ref);
            record -> ref = strdup(tok);
        } else if (i == 4) {
            for (int j = 0; j < strlen(tok); j++)
                if (tok[j] == ',')
                    numAlleles++;
            if (record -> alts != NULL)
                free(record -> alts);
            record -> alts = strdup(tok);
        } else if (i > 8) {
            record -> genotypes[i - 9].left = MISSING;
            record -> genotypes[i - 9].right = MISSING;
            record -> genotypes[i - 9].isPhased = false;
            char* start = strdup(tok);
            char* next = start;
            if (start[0] != '.')
                record -> genotypes[i - 9].left = (int) strtol(start, &next, 10);
            if (next[0] == '|')
                record -> genotypes[i - 9].isPhased = true;
            if ((next[0] == '|' || next[0] == '/') && next[1] != '.')
                record -> genotypes[i - 9].right = (int) strtol(next + 1, (char**) NULL, 10);
        }
    }
    record -> numAlleles = numAlleles;

    // Successfully parsed record.
    return true;

}

void parse_vcf(Replicate_t* replicate, InputStream_t* inputStream) {
    int recordIndex = 0;
    while (true) {

        // Allocate memory for a new record.
        Record_t* record = (Record_t*) calloc(1, sizeof(Record_t));
        record -> genotypes = (Genotype_t*) calloc(replicate -> numSamples, sizeof(Genotype_t));
        record -> numSamples = replicate -> numSamples;
        record -> recordIndex = recordIndex;

        // While not EOF, add record to the list.
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
        recordIndex++;
    }

}

Replicate_t* parse_ms(InputStream_t* inputStream, int length, bool hap) {
    if (inputStream == NULL)
        return NULL;

    // Eat lines until segsites encountered.
    do {
        ks_getuntil(inputStream -> fpIn, '\n', inputStream -> buffer, 0);
    } while (!ks_eof(inputStream -> fpIn) && strncmp(inputStream -> buffer -> s, "segsites:", 9) != 0);

    // If EOF, no replicate can be parsed. Return NULL.
    if (ks_eof(inputStream -> fpIn))
        return NULL;

    // List of the segregating site positions.
    kvec_t(int) positions;
    kv_init(positions);

    // List of the lineages in the ms replicate.
    kvec_t(char*) lineages;
    kv_init(lineages);

    // Get number of segsites.
    int segsites = (int) strtol(inputStream -> buffer -> s + 10, (char**) NULL, 10); 

    // Eat lines until positions encountered.
    do {
        ks_getuntil(inputStream -> fpIn, '\n', inputStream -> buffer, 0);
    } while (strncmp(inputStream -> buffer -> s, "positions:", 10) != 0);

    // Get the position of each segsite and add to list.
    int numSpaces = 0, prevIndex = 0, prevPos = -1;
    for (int i = 0; i <= inputStream -> buffer -> l; i++) {
        if (inputStream -> buffer -> s[i] == ' ' || i == inputStream -> buffer -> l) {
            if (numSpaces > 0) {
                // Convert position to base pairs. Ensure each position in bp is unique.
                int pos = (int) (length * strtod(inputStream -> buffer -> s + prevIndex, (char**) NULL));
                if (pos == prevPos)
                    pos += 1;
                prevPos = pos;
                kv_push(int, positions, pos); 
            }
            prevIndex = i;
            numSpaces++;
        }
    }

    // Read in the lineages.
    int numLineages = 0;
    while (ks_getuntil(inputStream -> fpIn, '\n', inputStream -> buffer, 0) > 0 && strncmp(inputStream -> buffer -> s, "//", 2) != 0) {
        kv_push(char*, lineages, strndup(inputStream -> buffer -> s, inputStream -> buffer -> l));
        numLineages++;
    }

    // Create our replicate. Each sample is diploid.
    Replicate_t* replicate = (Replicate_t*) calloc(1, sizeof(Replicate_t));
    replicate -> numSamples = hap ? numLineages : numLineages / 2;
    replicate -> numRecords = segsites;

    // For each segsite.
    for (int i = 0; i < replicate -> numRecords; i++) {
        // Create record.
        Record_t* temp = (Record_t*) calloc(1, sizeof(Record_t));
        temp -> genotypes = (Genotype_t*) calloc(replicate -> numSamples, sizeof(Genotype_t));
        temp -> position = kv_A(positions, i);
        temp -> numSamples = replicate -> numSamples; 
        temp -> numAlleles = 2;
        temp -> recordIndex = i;
        for (int j = 0; j < replicate -> numSamples; j++) {
            // MS lineages are phased.
            temp -> genotypes[j].isPhased = true;
            if (hap) {
                temp -> genotypes[j].left = kv_A(lineages, j)[i] - '0'; 
                temp -> genotypes[j].right = MISSING;
            } else {
                // Same pair of lineages creates a diploid sample.
                temp -> genotypes[j].left = kv_A(lineages, 2 * j)[i] - '0'; 
                temp -> genotypes[j].right = kv_A(lineages, 2 * j + 1)[i] - '0'; 
            }
        }
        // Add record to replicate list.
        if (replicate -> headRecord == NULL) {
            replicate -> headRecord = temp;
            replicate -> tailRecord = temp;
        } else {
            replicate -> tailRecord -> nextRecord = temp;
            replicate -> tailRecord = temp;
        }
    }

    // Destroy positions and lineages lists.
    kv_destroy(positions);
    for (int i = 0; i < kv_size(lineages); i++) {
        free(kv_A(lineages, i));
    }
    kv_destroy(lineages);

    return replicate;
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