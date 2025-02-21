
// File: main.c
// Date: 15 February 2025
// Version 1:
// Author: T. Quinn Smith
// Principal Investigator: Zachary A. Szpiech
// Purpose: Main functionality for EGGS.

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>
#include "../lib/ketopt.h"
#include "../lib/kvec.h"
#include "../lib/kstring.h"
#include "../lib/kseq.h"
#include "../lib/zlib.h"

// We use kseq to read in from stdin.
#define BUFFER_SIZE 4096
KSTREAM_INIT(gzFile, gzread, BUFFER_SIZE)

// We want random floats between [0, 1).
#define rand() ((float) rand() / (float) (RAND_MAX))

// Print our booleans.
#define PRINT_BOOL(X) (X ? "true" : "false")

// Operate on VCF input file.
// Accepts:
//  gzFile fp_out -> The file we are printing to.
//  All other arguments are the same from calling function.
// Returns: void.
void vcf(gzFile fp_out, kstream_t* fp_in, kstring_t* buffer, char* out, double missing, bool unphased, bool unpolarized) {
    // Echo the header information.
    while (!ks_eof(fp_in) && strncmp(ks_str(buffer), "#CHROM", 5) != 0) { 
        gzprintf(fp_out, "%s\n", ks_str(buffer)); 
        ks_getuntil(fp_in, '\n', buffer, 0); 
    } 
    if (ks_eof(fp_in)) return; 
    gzprintf(fp_out, "##eggs=<missing=%lf,unphased=%s,unpolarized=%s,out=%s>\n", missing, PRINT_BOOL(unphased), PRINT_BOOL(unpolarized), out); 
    gzprintf(fp_out, "%s\n", ks_str(buffer)); 
    kstring_t* temp = NULL; 
    kstring_t* leftGeno = init_kstring(NULL); 
    kstring_t* rightGeno = init_kstring(NULL); 
    int numAlts, slashIndex, colonIndex, numTok; 
    bool switchStates;
    while(true) { 
        ks_getuntil(fp_in, '\n', buffer, 0); 
        if (ks_eof(fp_in)) 
            break; 
        numAlts = 1; slashIndex = -1; colonIndex = 0, numTok = 1; 
        switchStates = false; 
        // VCF files are tab delimited.
        char* token = strtok(ks_str(buffer), "\t"); 
        while (token != NULL) { 
            // Count the number of alt alleles.
            if (numTok == 5) { 
                for (int j = 0; j < strlen(token); j++) 
                    if (token[j] == ',') 
                        numAlts++; 
                if (numAlts == 2 && unpolarized && rand() < 0.5) 
                    switchStates = true; 
                gzprintf(fp_out, "\t%s", token);
            // Parse the genotypes.
            } else if (numTok > 9) { 
                // Locate the '/' or '|' and ':' in the genotype.
                colonIndex = strlen(token); 
                for (int j = 0; j < strlen(token); j++) { 
                    if (token[j] == '/' || token[j] == '|')
                        slashIndex = j; 
                    if (token[j] == ':') { 
                        colonIndex = j; 
                        break; 
                    } 
                } 
                // Get the left and right allele.
                ks_overwriten(token, slashIndex, leftGeno); 
                ks_overwriten(token + slashIndex + 1, (colonIndex - slashIndex), rightGeno); 
                // Note, we do not switch the REF and ALT alleles to unpolarize because '.' is a valid ALT but 
                //  not a valid REF. For the same reason, we compare characters and do not use bit ops.
                if (switchStates) { 
                    if (ks_str(leftGeno)[0] == '0') 
                        ks_str(leftGeno)[0] = '1'; 
                    if (ks_str(rightGeno)[0] == '0') 
                        ks_str(rightGeno)[0] = '1'; 
                    if (ks_str(leftGeno)[0] == '1') 
                        ks_str(leftGeno)[0] = '0'; 
                    if (ks_str(rightGeno)[0] == '1') 
                        ks_str(rightGeno)[0] = '0'; 
                } 
                gzprintf(fp_out, "\t");
                // Same logic as ms.
                if (unphased) { 
                    if (rand() < 0.5) { temp = leftGeno; leftGeno = rightGeno; rightGeno = leftGeno; } 
                    if (missing > 0) { 
                        if (rand() < missing) gzprintf(fp_out, "."); else gzprintf(fp_out, "%s", ks_str(leftGeno)); 
                        if (rand() < missing) gzprintf(fp_out, "/."); else gzprintf(fp_out, "/%s", ks_str(rightGeno)); 
                    } else { gzprintf(fp_out, "%s/%s", ks_str(leftGeno), ks_str(rightGeno)); } 
                } else { 
                    if (missing > 0) { 
                        if (rand() < missing) gzprintf(fp_out, "."); else gzprintf(fp_out, "%s", ks_str(leftGeno)); 
                        if (rand() < missing) gzprintf(fp_out, "|."); else gzprintf(fp_out, "|%s", ks_str(rightGeno)); 
                    } else { gzprintf(fp_out, "%s|%s", ks_str(leftGeno), ks_str(rightGeno)); } 
                } 
                if (token[colonIndex] == ':') 
                    gzprintf(fp_out, "%s", token + colonIndex); 
            } else { 
                if (numTok == 1) 
                    gzprintf(fp_out, "%s", token); 
                else 
                    gzprintf(fp_out, "\t%s", token); 
            } 
            numTok++; 
            token = strtok(NULL, "\t");
        } 
        gzprintf(fp_out, "\n"); 
    } 
    destroy_kstring(leftGeno); 
    destroy_kstring(rightGeno); 
}

// Parse input in VCF format.
// Accepts:
//  kstream_t* fp_in -> The input.
//  kstring_t* buffer -> A buffer to read in the input.
//  double missing -> Probability of a missing genotype.
//  char* out -> The output basename.
//  bool unphased -> If sites should be unphased.
//  bool unpolarized -> If biallelic sites should be unpolarized.
// Returns: void.
void parse_vcf(kstream_t* fp_in, kstring_t* buffer, double missing, char* out, bool unphased, bool unpolarized) {
    // If we are printing to a file or stdout.
    if (out) {
        kstring_t* outName = init_kstring(out);
        kputs(".vcf.gz", outName);
        gzFile fp_out = gzopen(ks_str(outName), "w");
        vcf(fp_out, fp_in, buffer, out, missing, unphased, unpolarized);
        gzclose(fp_out);
        destroy_kstring(outName);
    } else {
        out = "stdout";
        int fd_out = fileno(stdout);
        gzFile fp_out = gzdopen(fd_out, "w");
        vcf(fp_out, fp_in, buffer, out, missing, unphased, unpolarized);
        gzclose(fp_out);
    }
}

// Parse ms replicate to VCF file with haploid samples.
// Accepts:
//  gzFile fp_out -> The file we are printing to.
//  All other arguments are the same from calling function.
// Returns: void.
void haploid_ms(gzFile fp_out, char* out, int length, bool unphased, double missing, bool unpolarized, bool single, int numSegsites, int numSamples, kvec_t(double)* positions, kvec_t(kstring_t*)* lineages) {
    // Used to parse replicate.
    int prevPosition = 0, pos;
    char leftGeno, rightGeno, temp;
    bool switchStates = false; // If we unpolarize switch the allelic states.

    gzprintf(fp_out, "##fileformat=VCFv4.2\n"); 
    gzprintf(fp_out, "##eggs=<numSamples=%d,missing=%lf,unphased=%s,unpolarized=%s,single=%s,out=%s>\n", numSamples, missing, PRINT_BOOL(unphased), PRINT_BOOL(unpolarized), PRINT_BOOL(single), out); 
    gzprintf(fp_out, "##contig=<ID=chr1,length=%d>\n", length); 
    gzprintf(fp_out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"); 
    for (int i = 0; i < numSamples; i++) { 
        gzprintf(fp_out, "\ts%d", i); 
    } 
    gzprintf(fp_out, "\n"); 
    for (int i = 0; i < numSegsites; i++) { 
        pos = (int) (kv_A(*positions, i) * length); 
        if (pos == prevPosition) { pos += 1; } 
        prevPosition = pos; 
        gzprintf(fp_out, "chr1\t%d\t.\tA\tT\t.\t.\t.\t.", pos); 
        if (unpolarized && rand() < 0.5) { 
            switchStates = true; 
        } 
        for (int j = 0; j < numSamples; j++) { 
            leftGeno = kv_A(*lineages, j) -> s[i]; 
            if (switchStates) { 
                leftGeno ^= 1; 
            } 
            if (missing > 0 && rand() < missing) { 
                gzprintf(fp_out, "\t./."); 
            } else { 
                gzprintf(fp_out, "\t%c/.", leftGeno); 
            } 
        } 
        gzprintf(fp_out, "\n"); 
        switchStates = false; 
    } 
}

// Parse ms replicate to VCF file with diploid samples.
// Accepts:
//  gzFile fp_out -> The file we are printing to.
//  All other arguments are the same from calling function.
// Returns: void.
void diploid_ms(gzFile fp_out, char* out, int length, bool unphased, double missing, bool unpolarized, bool single, int numSegsites, int numSamples, kvec_t(double)* positions, kvec_t(kstring_t*)* lineages) {
    // Used to parse replicate.
    int prevPosition = 0, pos;
    char leftGeno, rightGeno, temp;
    bool switchStates = false; // If we unpolarize switch the allelic states.

    gzprintf(fp_out, "##fileformat=VCFv4.2\n"); 
    gzprintf(fp_out, "##eggs=<numSamples=%d,missing=%lf,unphased=%s,unpolarized=%s,single=%s,out=%s>\n", numSamples, missing, PRINT_BOOL(unphased), PRINT_BOOL(unpolarized), PRINT_BOOL(single), out); 
    gzprintf(fp_out, "##contig=<ID=chr1,length=%d>\n", length); 
    gzprintf(fp_out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"); 
    for (int i = 0; i < numSamples / 2; i++) { 
        gzprintf(fp_out, "\ts%d", i); 
    } 
    gzprintf(fp_out, "\n"); 
    for (int i = 0; i < numSegsites; i++) { 
        pos = (int) (kv_A(*positions, i) * length); 
        if (pos == prevPosition) { pos += 1; } 
        prevPosition = pos; 
        gzprintf(fp_out, "chr1\t%d\t.\tA\tT\t.\t.\t.\t.", pos); 
        // If we unpolarize the record, set plag.
        if (unpolarized && rand() < 0.5) { 
            switchStates = true; 
        } 
        for (int j = 0; j < numSamples / 2; j++) { 
            leftGeno = kv_A(*lineages, 2 * j) -> s[i]; 
            rightGeno = kv_A(*lineages, 2 * j + 1) -> s[i]; 
            // Unpolarize.
            if (switchStates) { 
                leftGeno ^= 1; 
                rightGeno ^= 1; 
            } 
            if (unphased) { 
                if (rand() < 0.5) { temp = leftGeno; leftGeno = rightGeno; rightGeno = leftGeno; } 
                if (missing > 0) { 
                    if (rand() < missing) gzprintf(fp_out, "\t."); else gzprintf(fp_out, "\t%c", leftGeno); 
                    if (rand() < missing) gzprintf(fp_out, "/."); else gzprintf(fp_out, "/%c", rightGeno); 
                } else { gzprintf(fp_out, "\t%c/%c", leftGeno, rightGeno); } 
            } else { 
                if (missing > 0) { 
                    if (rand() < missing) gzprintf(fp_out, "\t."); else gzprintf(fp_out, "\t%c", leftGeno); 
                    if (rand() < missing) gzprintf(fp_out, "|."); else gzprintf(fp_out, "|%c", rightGeno); 
                } else { gzprintf(fp_out, "\t%c|%c", leftGeno, rightGeno); } 
            } 
        } 
        gzprintf(fp_out, "\n"); 
        switchStates = false; 
    } 
}

// Prints ms replicate to VCF file.
// Accepts:
//  bool isEOFandOneRep -> If set, reached eof.
//  int numRep -> The current replicate number.
//  char* out -> The base output file name.
//  int length -> The length of the segment in bp.
//  bool unphased -> If set, the resulting output should be unphased.
//  bool unpolarized -> If set, the records will be randomly unpolarized.
//  bool single -> If set, then each sample contains one unphased genotype.
//  int numSegsites -> The number of segregating sites in the replicate.
//  int numSamples -> The number of numSamples to get from the replicate.
//  kvec_t(double) positions -> The list of segregating site positions.
//  kvec_t(kstring_t*) lineages -> The list of simulated lineages.
// Returns: void.
void parse_ms_replicate(bool isEOF, int numRep, char* out, int length, bool unphased, double missing, bool unpolarized, bool single, int numSegsites, int numSamples, kvec_t(double)* positions, kvec_t(kstring_t*)* lineages) {
    
    // If there is one replicate and no output filename, we print to stdout.
    if (isEOF && numRep == 1 && !out) {
        out = "stdout";
        gzFile fp_out = gzdopen(fileno(stdout), "w");
        if (single) {
            haploid_ms(fp_out, out, length, unphased, missing, unpolarized, single, numSegsites, numSamples, positions, lineages);
        } else {
            diploid_ms(fp_out, out, length, unphased, missing, unpolarized, single, numSegsites, numSamples, positions, lineages);
        }
        gzclose(fp_out);
        return;
    }

    // Write to the file
    kstring_t* outName = init_kstring(NULL);
    // Create the output file name.
    if (out) {
        outName -> s = calloc(strlen(out) + ((int) log10(numRep + 1.0) + 1) + 12, sizeof(char));
        sprintf(ks_str(outName), "%s_rep%d.vcf.gz\0", out, numRep);
    } else {
        // Without an output name we use "rep#"
        outName -> s = calloc(((int) log10(numRep + 1.0) + 1) + 12, sizeof(char));
        sprintf(ks_str(outName), "rep%d.vcf.gz\0", numRep);
    } 
    gzFile fp_out = gzopen(ks_str(outName), "w");
    if (single) {
        haploid_ms(fp_out, out, length, unphased, missing, unpolarized, single, numSegsites, numSamples, positions, lineages);
    } else {
        diploid_ms(fp_out, out, length, unphased, missing, unpolarized, single, numSegsites, numSamples, positions, lineages);
    }
    gzclose(fp_out);
    
    destroy_kstring(outName);
}

// Parse input in ms-style format.
// Accepts:
//  kstream_t* fp_in -> The input stream.
//  kstring_t* buffer -> A buffer to read in the input.
//  int length -> The length in base pairs of the region.
//  int numSamples -> The number of numSamples to output.
//  double missing -> Probability of a missing genotype.
//  char* out -> The output basename.
//  bool unphased -> If sites should be unphased.
//  bool unpolarized -> If biallelic sites should be unpolarized.
//  bool single -> If a sample should be treated as a single lineage.
// Returns: void.
void parse_ms(kstream_t* fp_in, kstring_t* buffer, int length, int numSamples, double missing, char* out, bool unphased, bool unpolarized, bool single) {

    // Initalize memory used to read in a replicate.
    kvec_t(double) positions;
	kv_init(positions);
    kvec_t(kstring_t*) lineages;
    kv_init(lineages);

    // Eat lines until "segsites:" is encountered.
    do {
        ks_getuntil(fp_in, '\n', buffer, 0);
    } while (!ks_eof(fp_in) && strncmp(ks_str(buffer), "segsites:", 9) != 0);

    if (ks_eof(fp_in))
        return;

    int numReplicate = 1;

    // Parse all of the replicates.
    while (true) {
        
        int segsites = (int) strtol((ks_str(buffer)) + 10, (char**) NULL, 10); 

        // Eat lines until "positions:" is encountered.
        do {
            ks_getuntil(fp_in, '\n', buffer, 0);
        } while (strncmp(ks_str(buffer), "positions:", 10) != 0);
        
        // Get the positions of the segsites.
        int numSpaces = 0, prevIndex;
        for (int i = 0; i <= ks_len(buffer); i++) {
            if (ks_str(buffer)[i] == ' ' || i == ks_len(buffer)) {
                if (numSpaces > 0) {
                    double pos = strtod((ks_str(buffer)) + prevIndex, (char**) NULL);
                    if (numSpaces > kv_size(positions)) {
                        kv_push(double, positions, pos); 
                    } else {
                        kv_A(positions, numSpaces - 1) = pos;
                    }
                }
                prevIndex = i;
                numSpaces++;
            }
        }

        // Now, read in all of the lineages.
        int numLineages = 0;
        while (ks_getuntil(fp_in, '\n', buffer, 0) > 0 && strncmp(ks_str(buffer), "segsites:", 9) != 0) {
            if (numLineages >= kv_size(lineages)) {
                kstring_t* temp = calloc(1, sizeof(kstring_t));
                kputs(ks_str(buffer), temp);
                kv_push(kstring_t*, lineages, temp); 
            } else {
                ks_overwrite(ks_str(buffer), kv_A(lineages, numLineages));
            }
            numLineages++;
        }

        do {
            ks_getuntil(fp_in, '\n', buffer, 0);
        } while (!ks_eof(fp_in) && strncmp(ks_str(buffer), "segsites:", 9) != 0);

        // Check to make sure the number of lineages satisfies the number of samples.
        //  If not set, use all lineages.
        if (single && numSamples > numLineages) {
            printf("Requested %d lineages to sample. Only %d exist. Exiting!\n", numSamples, numLineages);
            break;
        } else if (numSamples > numLineages / 2) {
            printf("Requested %d diploids to sample. Only %d exist. Exiting!\n", numSamples, numLineages / 2);
            break;
        } else if (numSamples == -1) {
            // Use all avalibale lineages.
            parse_ms_replicate(ks_eof(fp_in), numReplicate, out, length, unphased, missing, unpolarized, single, segsites, numLineages, &positions, &lineages);
        } else {
            // Use subset of lineages.
            parse_ms_replicate(ks_eof(fp_in), numReplicate, out, length, unphased, missing, unpolarized, single, segsites, numSamples, &positions, &lineages);
        }

        numReplicate++;
        
        if (ks_eof(fp_in)) {
            break;
        }

    }

    // Free memory.
    kv_destroy(positions);
    for (int i = 0; i < kv_size(lineages); i++) {
        free(ks_str(kv_A(lineages, i)));
    }
    kv_destroy(lineages);

}

// Print the help menu for EGGS.
// Accepts: void.
// Returns: void.
void print_help() {
    printf("\n");
    printf("EGGS v1.0 February 2025\n");
    printf("----------------------\n\n");
    printf("Written by T. Quinn Smith\n");
    printf("Principal Investigator: Zachary A. Szpiech\n");
    printf("The Pennsylvania State University\n\n");
    printf("Usage: eggs [-upsch] [-l length] [-n number_of_numSamples] [-m prob_missing_allele] [-o out_basename]\n\n");
    printf("Options:\n");
    printf("    -l,--length     INT             Sets length of segment in number of base pairs for ms-replicates.\n");
    printf("                                        Default 1,000,000.\n");
    printf("    -n,--numSamples INT             Output the first n samples.\n");
    printf("                                        Default maximum possible per replicate.\n");
    printf("    -m,--missing    DOUBLE          Genotypes are missing with supplied probability.\n");
    printf("                                        Default 0.\n");
    printf("    -o,--out        STR             Basename for output files.\n");
    printf("                                        Default \"rep\" prefix for multiple ms-style replicates.\n");
    printf("    -u,--unphased                   Left and right genotypes are swapped with a probability of 0.5\n");
    printf("    -p,--unpolarized                Biallelic site alleles swapped with a probability of 0.5\n");
    printf("    -s,--single                     Each sample contains one lineage.\n");
    printf("                                        ms-style input only. Ignores -u.\n");
    printf("    -h,--help                       Print help.\n");
    printf("\n");
}

static ko_longopt_t long_options[] = {
    {"unphased",            ko_no_argument,         'u'},
    {"unpolarized",         ko_no_argument,         'p'},
    {"single",              ko_no_argument,         's'},
    {"stdout",              ko_no_argument,         'c'},
    {"length",              ko_required_argument,   'l'},
    {"numSamples",          ko_required_argument,   'n'},
    {"missing",             ko_required_argument,   'm'},
    {"out",                 ko_required_argument,   'o'},
    {"help",                ko_no_argument,         'h'},
    {NULL, 0, 0}
};

int main(int argc, char** argv) {

    // Seed random number generator.
    srand(time(NULL));

    // Single character aliases for long options.
    const char *opt_str = "huvspcl:n:m:o:";
    ketopt_t options = KETOPT_INIT;
    int c;

    // Options.
    int length = 1000000;
    int numSamples = -1;
    double missing = 0;
    char* out = NULL;
    bool unphased = false;
    bool unpolarized = false;
    bool single = false;

    // Parse options.
    while ((c = ketopt(&options, argc, argv, 1, opt_str, long_options)) >= 0) {
        if (c == 'l') length = atoi(options.arg);
        else if (c == 'u') unphased = true;
        else if (c == 'm') missing = atof(options.arg);
        else if (c == 'n') numSamples = atoi(options.arg);
        else if (c == 'p') unpolarized = true;
        else if (c == 's') single = true;
        else if (c == 'o') out = strdup(options.arg);
        else if (c == 'h') { print_help(); return 0; }
        else if (c ==':') { fprintf(stderr, "Error! Option %s is missing an argument! Exiting ...\n", argv[options.i - 1]); return 1; }
        else { fprintf(stderr, "Error! \"%s\" is unknown! Exiting ...\n", argv[options.i - 1]); return 1; }
	}

    // Check valid missing probability.
    if (missing < 0) {
        fprintf(stderr, "Missing probability must be >= 0! Exiting ...\n");
        if (out) free(out);
        return 1;
    }

    // Check valid region length.
    if (length < 1000) {
        fprintf(stderr, "Length for region should be >= 1000! Exiting ...\n");
        if (out) free(out);
        return 1;
    }

    // We read from stdin.
    int f_in = fileno(stdin);
    gzFile file = gzdopen(f_in, "r");
    kstream_t* fp_in = ks_init(file);
    kstring_t* buffer = init_kstring(NULL);

    // If the first line signals a VCF file, then we parse it as a VCF.
    //  Otherwise, we assume it is ms-style.
    ks_getuntil(fp_in, '\n', buffer, 0);
    if (strncmp(ks_str(buffer), "##fileformat=VCF", 16) == 0) {
        parse_vcf(fp_in, buffer, missing, out, unphased, unpolarized);
    } else {
        parse_ms(fp_in, buffer, length, numSamples, missing, out, unphased, unpolarized, single);
    }

    // Free memory.
    destroy_kstring(buffer);
    ks_destroy(fp_in);
    gzclose(file);
    // Free output base name if used.
    if (out) free(out);
}