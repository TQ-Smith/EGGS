
// File: Main.c
// Date: 13 May 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Main logic of EGGS.

#include "Interface.h"
#include "GenotypeParser.h"
#include "kstring.h"
#include "Missingness.h"
#include <math.h>
#include <time.h>
#include "compact_bitset.h"

// Swap integer values.
#define SWAP(a, b, temp) { a = temp; a = b; b = temp; }

// Easy random uniform float with [0, 1)
#define rand() ((float) rand() / (float) (RAND_MAX))

// Print a simple VCF header from a replicate.
void print_vcf_header(Replicate_t* replicate, EggsConfig_t* eggsConfig, gzFile fpOut) {
    gzprintf(fpOut, "##fileformat=VCFv4.2\n");
    gzprintf(fpOut, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"); 
    // If -r was used, explicitly print values.
    if (eggsConfig -> betaMissing != NULL)
        gzprintf(fpOut, "##-b=%lf,%lf\n", eggsConfig -> meanMissing, eggsConfig -> stdMissing);
    // Print user command.
    gzprintf(fpOut, "##eggsCommand=%s\n", eggsConfig -> command);
    gzprintf(fpOut, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"); 
    // Print sample names if exits.
    if (replicate -> sampleNames == NULL)
        for (int i = 0; i < replicate -> numSamples; i++)
            gzprintf(fpOut, "\ts%d", i); 
    // Otherwise, an ms-replicate was provided and we generate simple sample names by index.
    else
        for (int i = 0; i < replicate -> numSamples; i++)
            gzprintf(fpOut, "\t%s", replicate -> sampleNames[i]);
    gzprintf(fpOut, "\n");
}


// Print a simple ms header from a replicate.
void print_ms_header(Replicate_t* replicate, EggsConfig_t* eggsConfig, gzFile fpOut) {
    gzprintf(fpOut, "eggsv1.0 \n\n");
}



// Prints a record.
// Accepts:
//  Record_t* record -> The record we are printing out.
//  Mask_t* mask -> If index contains MISSING, then the sample's genotype is masked out.
//  EggsConfig_t* eggsConfig -> Defines the logic of the record manipulation.
//  gzFile fpOut -> The output stream.
// Returns: void.
void print_record(Record_t* record, Mask_t* mask, EggsConfig_t* eggsConfig, gzFile fpOut) {
    // If ms-replicate, then we just use chrom 1.
    if (record -> chrom == NULL) 
        gzprintf(fpOut, "chr1");
    else 
        gzprintf(fpOut, "%s", record -> chrom);

    // Print position.
    gzprintf(fpOut, "\t%d\t.", record -> position);

    // If ms-replicate, assign REF a dummy allele.
    if (record -> ref == NULL)
        gzprintf(fpOut, "\tA");
    else 
        gzprintf(fpOut, "\t%s", record -> ref);

    // If ms-replicate, assign ALTS a dummy allele.
    if (record -> alts == NULL)
        gzprintf(fpOut, "\tT");
    else 
        gzprintf(fpOut, "\t%s", record -> alts);

    // Empty fields.
    gzprintf(fpOut, "\t.\t.\t.\tGT");

    // Now, we print the genotypes.

    // Used to swap allele to unphase.
    int tempInt = 0;
    // Flag is set if biallelic site should be unpolarized with 50/50 chance.
    bool swapStates = false;
    if (eggsConfig -> unpolarize && record -> numAlleles == 2 && rand() < 0.5)
        swapStates = true;

    // Flag to set if biallelic site should be deanimated.
    bool isTransition = false;
    if (eggsConfig -> probTransition != 0 && record -> numAlleles == 2 && rand() < eggsConfig -> probTransition)
        isTransition = true;

    for (int i = 0; i < record -> numSamples; i++) {
        // If the genotype is masked out, then skip rest of logic.
        if (mask != NULL && mask -> missing[record -> recordIndex][i] == MISSING) {
            record -> genotypes[i].isPhased = false;
            record -> genotypes[i].left = MISSING;
            record -> genotypes[i].right = MISSING;
        } else {

            // If genotypes should be unphased with a 50/50 chance.
            if (eggsConfig -> unphase && !eggsConfig -> pseudohap) {
                record -> genotypes[i].isPhased = false;
                if (rand() < 0.5)
                    SWAP(record -> genotypes[i].left, record -> genotypes[i].right, tempInt);
            }

            // If biallelic and locus should be unpolarized.
            if (swapStates && record -> genotypes[i].left != MISSING)
                record -> genotypes[i].left ^= 1;
            if (swapStates && record -> genotypes[i].right != MISSING)
                record -> genotypes[i].right ^= 1;

            // If the site should be demainated.
            if (isTransition && record -> genotypes[i].left == 0 && rand() < eggsConfig -> probDeamination)
                record -> genotypes[i].left = 1;
            if (isTransition && record -> genotypes[i].right == 0 && rand() < eggsConfig -> probDeamination)
                record -> genotypes[i].right = 1;

            // If pseudohap is set, then randomly pick one allele.
            if (eggsConfig -> pseudohap) {
                record -> genotypes[i].isPhased = false;
                if (rand() < 0.5)
                    record -> genotypes[i].right = record -> genotypes[i].left;
                else
                    record -> genotypes[i].left = record -> genotypes[i].right;
            }

        }

        // Print out sample's genotype.
        if (record -> genotypes[i].left == MISSING)
            gzprintf(fpOut, "\t.");
        else 
            gzprintf(fpOut, "\t%d", record -> genotypes[i].left);
        if (record -> genotypes[i].isPhased)
            gzprintf(fpOut, "|");
        else 
            gzprintf(fpOut, "/");
        if (record -> genotypes[i].right == MISSING)
            gzprintf(fpOut, ".");
        else 
            gzprintf(fpOut, "%d", record -> genotypes[i].right);
    }
    gzprintf(fpOut, "\n");
}



// Print the list of records in a replicate, assuming output is in ms format.
// Accepts:
//  Replicate_t* replicate -> The records to print.
//  Mask_t* mask -> The mask to apply.
//  EggsConfig_t* eggsConfig -> The user parameters.
//  gzFile fpOut -> The output stream.
// Returns: void.
void print_replicate_ms(Replicate_t* replicate, Mask_t* mask, EggsConfig_t* eggsConfig, gzFile fpOut) {

    int n = replicate -> numSamples * 2;         // number of bitsets
    int bitlen = replicate -> numRecords;  // each bitset has bitlen bits
    
    gzprintf(fpOut, "//\n");
    gzprintf(fpOut, "segsites: %d\n", replicate -> numRecords);
    gzprintf(fpOut, "positions: ");
    
    Record_t* record = replicate -> headRecord;
    for (int i = 0; i < replicate -> numRecords; i++) {
        if(record -> position > eggsConfig->length){
            fprintf(stderr, "ERROR: positions must be less than or equal to the segment length provided with -l.\n");
            fprintf(stderr, "Found position %d, but segment length is %d.\n", record -> position, eggsConfig -> length);
            exit(EXIT_FAILURE);
        }
        gzprintf(fpOut, "%.8f ", record -> position*1.0 / eggsConfig -> length);
        record = record -> nextRecord; //next record
    }

    gzprintf(fpOut, "\n");

    CompactBitset* bitsets[n];
    for (int i = 0; i < n; ++i) {
        bitsets[i] = cb_create(bitlen);
    }
    
    record = replicate -> headRecord;
    for (int p = 0; p < replicate -> numRecords; p++) {
        // Used to swap allele to unphase.
        int tempInt = 0;
        // Flag is set if biallelic site should be unpolarized with 50/50 chance.
        bool swapStates = false;
        if (eggsConfig -> unpolarize && record -> numAlleles == 2 && rand() < 0.5){
            swapStates = true;
        }

        // Flag to set if biallelic site should be deanimated.
        bool isTransition = false;
        if (eggsConfig -> probTransition != 0 && record -> numAlleles == 2 && rand() < eggsConfig -> probTransition)
            isTransition = true;

        // do the processing in file
        for (int i = 0; i < record -> numSamples; i++) {
            // If the genotype is masked out, then skip rest of logic.
            if (mask != NULL && mask -> missing[record -> recordIndex][i] == MISSING) {
                fprintf(stderr, "ERROR: missing entries not allowed in ms output\n");
                exit(EXIT_FAILURE);  
            } else {
                // If genotypes should be unphased with a 50/50 chance.
                if (eggsConfig -> unphase && !eggsConfig -> pseudohap) {
                    record -> genotypes[i].isPhased = false;
                    if (rand() < 0.5)
                        SWAP(record -> genotypes[i].left, record -> genotypes[i].right, tempInt);
                }

                // If biallelic and locus should be unpolarized.
                if (swapStates && record -> genotypes[i].left != MISSING)
                    record -> genotypes[i].left ^= 1;
                if (swapStates && record -> genotypes[i].right != MISSING)
                    record -> genotypes[i].right ^= 1;

                // If the site should be demainated.
                if (isTransition && record -> genotypes[i].left == 0 && rand() < eggsConfig -> probDeamination)
                    record -> genotypes[i].left = 1;
                if (isTransition && record -> genotypes[i].right == 0 && rand() < eggsConfig -> probDeamination)
                    record -> genotypes[i].right = 1;

                // If pseudohap is set, then randomly pick one allele.
                if (eggsConfig -> pseudohap) {
                    record -> genotypes[i].isPhased = false;
                    if (rand() < 0.5)
                        record -> genotypes[i].right = record -> genotypes[i].left;
                    else
                        record -> genotypes[i].left = record -> genotypes[i].right;
                }
            }

            // Print out sample's genotype.
            if (record -> genotypes[i].left == MISSING || record -> genotypes[i].right == MISSING) {
                fprintf(stderr, "ERROR: missing entries not allowed in ms output\n");
                exit(EXIT_FAILURE);  
            }

            if(record -> genotypes[i].left == 1 )
                cb_set_bit(bitsets[2*i], p);

            if(record -> genotypes[i].right == 1 )
                cb_set_bit(bitsets[2*i+1], p);
        }
        record = record -> nextRecord; //next record
    }
        
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < bitsets[i]->nbits; ++k) {
            gzprintf(fpOut, "%d", cb_get_bit(bitsets[i], k));
        }
        gzprintf(fpOut, "\n");
    }

    // Cleanup
    for (int i = 0; i < n; ++i) {
        cb_destroy(bitsets[i]);
    }
}


// Print the list of records in a replicate.
// Accepts:
//  Replicate_t* replicate -> The records to print.
//  Mask_t* mask -> The mask to apply.
//  EggsConfig_t* eggsConfig -> The user parameters.
//  gzFile fpOut -> The output stream.
// Returns: void.
void print_replicate(Replicate_t* replicate, Mask_t* mask, EggsConfig_t* eggsConfig, gzFile fpOut) {
    Record_t* temp = replicate -> headRecord;
    for (int i = 0; i < replicate -> numRecords; i++) {
        if (mask == NULL)
            print_record(temp, NULL, eggsConfig, fpOut);
        else 
            // Pass the record's corresponding mask.
            print_record(temp, mask, eggsConfig, fpOut);
        temp = temp -> nextRecord;
    }
}

// Get mu and sigma of missing genotypes per site.
// Accepts:
//  InputStream_t* inputStream -> The VCF file to read.
//  double* mu -> Sets the mean per site.
//  double* sigma -> Sets the sigma per site.
// Returns: void.
void get_mu_sigma(InputStream_t* inputStream, double* mu, double* sigma) {

    // Read the VCF file.
    Replicate_t* replicate = init_vcf_replicate(inputStream);
    parse_vcf(replicate, inputStream);

    MissingDistribution_t* dis = init_missing_distribution(replicate);
    
    // Calculate and set mean.
    double mean = 0;
    for (int i = 0; i < dis -> numRecords; i++)
        mean += dis -> proportions[i] / (double) dis -> numRecords;
    *mu = mean;

    // Calculate and set variance.
    double var = 0;
    for (int i = 0; i < dis -> numRecords; i++)
        var += (dis -> proportions[i] - mean) * (dis -> proportions[i] - mean) / (dis -> numRecords - 1.0);
    *sigma = sqrt(var);
    
    destroy_missing_distribution(dis);
    destroy_replicate(replicate);
}

int main(int argc, char* argv[]) {

    srand(time(NULL));

    // Get CLI configuration.
    EggsConfig_t* eggsConfig = init_eggs_configuration(argc, argv);
    if (eggsConfig == NULL)
        return -1;

    // If a VCF file was given for random missingness, calculate mean and standard deviation per site.
    if (eggsConfig -> betaMissing != NULL && eggsConfig -> meanMissing == -1 && eggsConfig -> stdMissing == -1) {
        InputStream_t* inputStream = init_input_stream(eggsConfig -> betaMissing);
        double mu = 0, sigma = 0;
        get_mu_sigma(inputStream, &mu, &sigma);
        eggsConfig -> meanMissing = mu;
        eggsConfig -> stdMissing = sigma;
        destroy_input_stream(inputStream);
        if (eggsConfig -> meanMissing >= 1 || eggsConfig -> meanMissing <= 0 || eggsConfig -> stdMissing <= 0 || eggsConfig -> stdMissing * eggsConfig -> stdMissing >= eggsConfig -> meanMissing * (1 - eggsConfig -> meanMissing)) {
            fprintf(stderr, "mu=%lf,sigma=%lf must satisfy parameters for a beta distribution. Exiting!\n", eggsConfig -> meanMissing, eggsConfig -> stdMissing);
            destroy_eggs_configuration(eggsConfig);
            return -1;
        }
    }

    // Read from stdin.
    InputStream_t* inputStream = init_input_stream(NULL);

    // Replicate counter for ms-style replicates.
    int numReps = 1;

    // Read first line from stdin to determine if it is a VCF file or ms-style input.
    ks_getuntil(inputStream -> fpIn, '\n', inputStream -> buffer, 0);

    // If there is not stdin to read from, then we exit.
    if (ks_eof(inputStream -> fpIn)) {
        destroy_input_stream(inputStream);
        destroy_eggs_configuration(eggsConfig);
        return -1;
    }
    
    // If VCF.
    if (strncmp(inputStream -> buffer -> s, "##fileformat=VCF", 16) == 0) {
        // If output basename was given, open file. Otherwise, we are printing to stdout.
        gzFile fpOut;
        if (eggsConfig -> outFile != NULL) {
            kstring_t* outName = (kstring_t*) calloc(1, sizeof(kstring_t));
            ksprintf(outName, "%s.vcf.gz", eggsConfig -> outFile);
            fpOut = gzopen(outName -> s, "w");
            free(outName -> s); free(outName);
        } else {
            fpOut = gzdopen(fileno(stdout), "w");
        }

        Mask_t* mask = NULL;
        Replicate_t* replicate = init_vcf_replicate(inputStream);
        parse_vcf(replicate, inputStream);
        
        if (eggsConfig -> maskFile != NULL || eggsConfig -> randomMissing != NULL) {

            // Read in mask file.
            InputStream_t* maskInput = init_input_stream(eggsConfig -> maskFile != NULL ? eggsConfig -> maskFile : eggsConfig -> randomMissing);
            Replicate_t* maskReplicate = init_vcf_replicate(maskInput);
            parse_vcf(maskReplicate, maskInput);

            // Calculate proportion of missing samples at each site.
            MissingDistribution_t* dis = init_missing_distribution(maskReplicate);
            if (eggsConfig -> maskFile != NULL)
                mask = create_missing_mask(dis, replicate -> numSamples, replicate -> numRecords);
            else 
                mask = create_random_mask(dis, replicate -> numSamples, replicate -> numRecords, -1, -1);
            
            destroy_input_stream(maskInput);
            destroy_replicate(maskReplicate);
            destroy_missing_distribution(dis);
        // If we want to create a beta mask.
        } else if (eggsConfig -> betaMissing != NULL) {
            // Create our beta mask.
            mask = create_random_mask(NULL, replicate -> numSamples, replicate -> numRecords, eggsConfig -> meanMissing, eggsConfig -> stdMissing);

        } 

        if(eggsConfig -> ms) {
            print_ms_header(replicate, eggsConfig, fpOut);
        }else{
            print_vcf_header(replicate, eggsConfig, fpOut);
        }

        // If no mask was created, then we can just sequentually print records out.
        if (mask == NULL) {
            Record_t* record = (Record_t*) calloc(1, sizeof(Record_t));
            record -> genotypes = (Genotype_t*) calloc(replicate -> numSamples, sizeof(Genotype_t));
            record -> numSamples = replicate -> numSamples;
            while (get_next_vcf_record(record, inputStream))
                print_record(record, NULL, eggsConfig, fpOut);
            destroy_record(record);
        } else {
            print_replicate(replicate, mask, eggsConfig, fpOut);
        }
        
        destroy_replicate(replicate);
        gzclose(fpOut);
        destroy_mask(mask);

    // Otherwise, parse as ms-style input.
    } else {
        Replicate_t* replicate = NULL;
        MissingDistribution_t* dis = NULL;

        if (eggsConfig -> maskFile != NULL || eggsConfig -> randomMissing != NULL) {
            InputStream_t* maskInput = init_input_stream(eggsConfig -> maskFile != NULL ? eggsConfig -> maskFile : eggsConfig -> randomMissing);
            Replicate_t* maskReplicate = init_vcf_replicate(maskInput);
            parse_vcf(maskReplicate, maskInput);
            dis = init_missing_distribution(maskReplicate);
            destroy_input_stream(maskInput);
            destroy_replicate(maskReplicate);
        }



        // For each replicate in the ms-style input.
        while ((replicate = parse_ms(inputStream, eggsConfig -> length, eggsConfig -> hap)) != NULL) {
            Mask_t* mask = NULL;
            if (eggsConfig -> maskFile != NULL)
                mask = create_missing_mask(dis, replicate -> numSamples, replicate -> numRecords);
            // Generate our beta mask.
            else if (eggsConfig -> betaMissing != NULL) 
                mask = create_random_mask(NULL, replicate -> numSamples, replicate -> numRecords, eggsConfig -> meanMissing, eggsConfig -> stdMissing);
            // If user defined distribution was given.
            else if (eggsConfig -> randomMissing != NULL) 
                mask = create_random_mask(dis, replicate -> numSamples, replicate -> numRecords, -1, -1);
            
            gzFile fpOut;
            // If no basename given and only one replicate, then we print to stdout.
            if (ks_eof(inputStream -> fpIn) && eggsConfig -> outFile == NULL && numReps == 1) {
                fpOut = gzdopen(fileno(stdout), "w");
            // If output basename was given, then print replicates to basename.
            } else if (eggsConfig -> outFile != NULL) {
                kstring_t* outName = (kstring_t*) calloc(1, sizeof(kstring_t));
                if(eggsConfig -> ms) {
                    // If ms-style output, then we print to basename with replicate number.
                    ksprintf(outName, "%s_rep%d.ms.gz", eggsConfig -> outFile, numReps);
                } else{
                    // If VCF output, then we print to basename with replicate number.
                    ksprintf(outName, "%s_rep%d.vcf.gz", eggsConfig -> outFile, numReps);
                }
                fpOut = gzopen(outName -> s, "w");
                free(outName -> s); free(outName);
            } else {
                // If no basename was given use "rep".
                kstring_t* outName = (kstring_t*) calloc(1, sizeof(kstring_t));

                if(eggsConfig -> ms) {
                    // If ms-style output, then we print to "rep" with replicate number.
                    ksprintf(outName, "rep%d.ms.gz", numReps);
                } else{
                    // If VCF output, then we print to "rep" with replicate number.
                    ksprintf(outName, "rep%d.vcf.gz", numReps);
                }
                fpOut = gzopen(outName -> s, "w");
                free(outName -> s); free(outName);
            }

            // Print out replicate.
            if(eggsConfig -> ms) {
                print_ms_header(replicate, eggsConfig, fpOut);
                print_replicate_ms(replicate, mask, eggsConfig, fpOut);
            }else{
                print_vcf_header(replicate, eggsConfig, fpOut);
                print_replicate(replicate, mask, eggsConfig, fpOut);
            }
            

            gzclose(fpOut);
            destroy_replicate(replicate);
            destroy_mask(mask);
            numReps++;
        }
        destroy_missing_distribution(dis);
    }

    destroy_eggs_configuration(eggsConfig);
    destroy_input_stream(inputStream);

    return 0;
}