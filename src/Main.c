
// File: Main.c
// Date: 13 May 2025
// Authors: T. Quinn Smith and Dr. Amatur Rahman
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Main logic of EGGS.

#include "Interface.h"
#include "GenotypeParser.h"
#include "kstring.h"
#include "Missingness.h"
#include <math.h>
#include <time.h>

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

// Print the replicate to ms format.
void print_ms_replicate(Replicate_t* replicate, EggsConfig_t* eggsConfig, gzFile fpOut) {
    gzprintf(fpOut, "\n//\n");
    gzprintf(fpOut, "segsites: %d\n", replicate -> numRecords);
    gzprintf(fpOut, "positions: ");

    // Get closest multiple of a 1000 if not default.
    int sizeOfSegment;
    if (eggsConfig -> length > 0) 
        sizeOfSegment = eggsConfig -> length;
    else
        sizeOfSegment = ((int) (replicate -> tailRecord -> position / 1000) + 1) * 1000;
    gzprintf(fpOut, "%.8lf", replicate -> headRecord -> position / (double) sizeOfSegment);
    for (Record_t* temp = replicate -> headRecord -> nextRecord; temp != NULL; temp = temp -> nextRecord)
        gzprintf(fpOut, " %.8lf", temp -> position / (double) sizeOfSegment);
    gzprintf(fpOut, "\n");

    // We process the sites first.
    for (Record_t* temp = replicate -> headRecord; temp != NULL; temp = temp -> nextRecord) {
        // Used to swap allele to unphase.
        int tempInt = 0;
        // Flag is set if biallelic site should be unpolarized with 50/50 chance.
        bool swapStates = false;
        if (eggsConfig -> unpolarize && rand() < 0.5)
            swapStates = true;

        // Flag to set if biallelic site should be deanimated.
        bool isTransition = false;
        if (eggsConfig -> probTransition != 0 && rand() < eggsConfig -> probTransition)
            isTransition = true;

        for (int i = 0; i < replicate -> numSamples; i++) {

            // Put non missing in left. If missing, use ancestral. If multiallelic use alternative.
            if (!temp -> genotypes[i].isPhased) {
                if (temp -> genotypes[i].left == MISSING && temp -> genotypes[i].right != MISSING)
                    SWAP(temp -> genotypes[i].left, temp -> genotypes[i].right, tempInt);
                if (temp -> genotypes[i].right == MISSING)
                    temp -> genotypes[i].right = 0;
            }
            if (temp -> genotypes[i].left > 0)
                temp -> genotypes[i].left = 1;
            if (temp -> genotypes[i].right > 0)
                temp -> genotypes[i].right = 1;

            // If genotypes should be unphased with a 50/50 chance.
            if (eggsConfig -> unphase && !eggsConfig -> pseudohap) {
                temp -> genotypes[i].isPhased = false;
                if (rand() < 0.5)
                    SWAP(temp -> genotypes[i].left, temp -> genotypes[i].right, tempInt);
            }

            if (swapStates && temp -> genotypes[i].left != MISSING) 
                temp -> genotypes[i].left ^= 1;
            if (swapStates && temp -> genotypes[i].right != MISSING)
                temp -> genotypes[i].right ^= 1;

            // If the site should be demainated.
            if (isTransition && temp -> genotypes[i].left == 0 && rand() < eggsConfig -> probDeamination)
                temp -> genotypes[i].left = 1;
            if (isTransition && temp -> genotypes[i].right == 0 && rand() < eggsConfig -> probDeamination)
                temp -> genotypes[i].right = 1;

            // If pseudohap is set, then randomly pick one allele.
            if (eggsConfig -> pseudohap) {
                temp -> genotypes[i].isPhased = false;
                if (rand() < 0.5)
                    temp -> genotypes[i].right = temp -> genotypes[i].left;
                else
                    temp -> genotypes[i].left = temp -> genotypes[i].right;
            }
        }
    }

    // Then we print the lineages.
    for (int i = 0; i < replicate -> numSamples; i++) {
        for (Record_t* temp = replicate -> headRecord; temp != NULL; temp = temp -> nextRecord)
            gzprintf(fpOut, "%d", temp -> genotypes[i].left);
        if (!eggsConfig -> hap)
            for (Record_t* temp = replicate -> headRecord; temp != NULL; temp = temp -> nextRecord)
                gzprintf(fpOut, "%d", temp -> genotypes[i].right);
        gzprintf(fpOut, "\n");
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
    if (!eggsConfig -> msOutput) {
        Record_t* temp = replicate -> headRecord;
        for (int i = 0; i < replicate -> numRecords; i++) {
            if (mask == NULL)
                print_record(temp, NULL, eggsConfig, fpOut);
            else 
                // Pass the record's corresponding mask.
                print_record(temp, mask, eggsConfig, fpOut);
            temp = temp -> nextRecord;
        }
    } else {
        print_ms_replicate(replicate, eggsConfig, fpOut);
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
            if (eggsConfig -> msOutput)
                ksprintf(outName, "%s.ms.gz", eggsConfig -> outFile);
            else
                ksprintf(outName, "%s.vcf.gz", eggsConfig -> outFile);
            fpOut = gzopen(outName -> s, "w");
            free(outName -> s); free(outName);
        } else {
            fpOut = gzdopen(fileno(stdout), "w");
        }

        Mask_t* mask = NULL;
        Replicate_t* replicate = init_vcf_replicate(inputStream);
        
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

        // If no ms output and mask was created, then we can just sequentually print records out.
        if (!eggsConfig -> msOutput && mask == NULL) {
            Record_t* record = (Record_t*) calloc(1, sizeof(Record_t));
            record -> genotypes = (Genotype_t*) calloc(replicate -> numSamples, sizeof(Genotype_t));
            record -> numSamples = replicate -> numSamples;
            while (get_next_vcf_record(record, inputStream))
                print_record(record, NULL, eggsConfig, fpOut);
            destroy_record(record);
        } else {
            parse_vcf(replicate, inputStream);
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
        int length;
        if (eggsConfig -> length <= 0 )
            length = 1000000;
        else 
            length = eggsConfig -> length;
        while ((replicate = parse_ms(inputStream, length, eggsConfig -> hap)) != NULL) {
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
                if (eggsConfig -> msOutput)
                    ksprintf(outName, "%s_rep%d.ms.gz", eggsConfig -> outFile, numReps);
                else
                    ksprintf(outName, "%s_rep%d.vcf.gz", eggsConfig -> outFile, numReps);
                fpOut = gzopen(outName -> s, "w");
                free(outName -> s); free(outName);
            // If no basename was given use "rep".
            } else {
                kstring_t* outName = (kstring_t*) calloc(1, sizeof(kstring_t));
                if (eggsConfig -> msOutput)
                    ksprintf(outName, "%s_rep%d.ms.gz", eggsConfig -> outFile, numReps);
                else
                    ksprintf(outName, "rep%d.vcf.gz", numReps);
                fpOut = gzopen(outName -> s, "w");
                free(outName -> s); free(outName);
            }

            // Print out replicate.
            if (eggsConfig -> msOutput)
                gzprintf(fpOut, "%s\n", eggsConfig -> command);
            else
                print_vcf_header(replicate, eggsConfig, fpOut);
            print_replicate(replicate, mask, eggsConfig, fpOut);

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