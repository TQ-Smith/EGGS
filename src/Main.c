
// File: Main.c
// Date: 13 May 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Main logic of EGGS.

#include "Interface.h"
#include "GenotypeParser.h"
#include "time.h"
#include "kstring.h"
#include "Missingness.h"
#include "kvec.h"

// Swap integer values.
#define SWAP(a, b, temp) { a = temp; a = b; b = temp; }

// Easy random uniform float with [0, 1)
#define rand() ((float) rand() / (float) (RAND_MAX))

// Print a simple VCF header from a replicate.
void print_vcf_header(Replicate_t* replicate, EggsConfig_t* eggsConfig, gzFile fpOut) {
    gzprintf(fpOut, "##fileformat=VCFv4.2\n"); 
    // If -r was used, explicitly print values.
    if (eggsConfig -> randomMissing != NULL)
        gzprintf(fpOut, "##-r=%lf,%lf\n", eggsConfig -> meanMissing, eggsConfig -> stdMissing);
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
    gzprintf(fpOut, "\t.\t.\t.\t.");

    // Now, we print the genotypes.

    // Used to swap allele to unphase.
    int tempInt = 0;
    // Flag is set if biallelic site should be unpolarized with 50/50 chance.
    bool swapStates = false;
    if (eggsConfig -> unpolarize && record -> numAlleles == 2 && rand() < 0.5)
        swapStates = true;

    for (int i = 0; i < record -> numSamples; i++) {
        // If the genotype is masked out, then skip rest of logic.
        if (mask != NULL && mask -> missing[i][record -> recordIndex] == MISSING) {
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
    Replicate_t* replicate = init_vcf_replicate(inputStream);
    Record_t* record = (Record_t*) calloc(1, sizeof(Record_t));
    record -> genotypes = (Genotype_t*) calloc(replicate -> numSamples, sizeof(Genotype_t));
    record -> numSamples = replicate -> numSamples;
            
    kvec_t(double) proportions;
    kv_init(proportions);

    double mean = 0;
    int numRecords = 0;

    // Calculate the proportion of missing samples for each locus.
    while (get_next_vcf_record(record, inputStream)) {
        int numMissing = 0;
        for (int i = 0; i < record -> numSamples; i++) 
            if (record -> genotypes[i].left == MISSING && record -> genotypes[i].right == MISSING)
                numMissing++;
        mean += numMissing / (double) record -> numSamples;
        kv_push(double, proportions, numMissing / (double) record -> numSamples); 
        numRecords++;
    }
    
    // Set mean and standard deviation.
    mean /= numRecords;
    *mu = mean;
    *sigma = 0;
    for (int i = 0; i < numRecords; i++)
        *sigma += (kv_A(proportions, i) - mean) * (kv_A(proportions, i) - mean) / (numRecords - 1);
    *sigma = sqrt(*sigma);

    destroy_record(record);
    destroy_replicate(replicate);
    kv_destroy(proportions);
}

int main(int argc, char* argv[]) {

    srand(time(NULL));

    // Get CLI configuration.
    EggsConfig_t* eggsConfig = init_eggs_configuration(argc, argv);
    if (eggsConfig == NULL)
        return -1;

    // If a VCF file was given for random missingness, calculate mean and standard deviation per site.
    if (eggsConfig -> randomMissing != NULL && eggsConfig -> meanMissing == -1 && eggsConfig -> stdMissing == -1) {
        InputStream_t* inputStream = init_input_stream(eggsConfig -> randomMissing);
        double mu = 0, sigma = 0;
        get_mu_sigma(inputStream, &mu, &sigma);
        eggsConfig -> meanMissing = mu;
        eggsConfig -> stdMissing = sigma;
        destroy_input_stream(inputStream);
        if (eggsConfig -> meanMissing >= 1 || eggsConfig -> meanMissing <= 0 || eggsConfig -> stdMissing <= 0 || eggsConfig -> stdMissing * eggsConfig -> stdMissing >= eggsConfig -> meanMissing * (1 - eggsConfig -> meanMissing)) {
            fprintf(stderr, "VCF mu=%lf,sigma=%lf must satisfy parameters for a beta distribution. Exiting!\n", eggsConfig -> meanMissing, eggsConfig -> stdMissing);
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

        // Read VCF header from stdin.
        Replicate_t* replicate = init_vcf_replicate(inputStream);
        print_vcf_header(replicate, eggsConfig, fpOut);

        // Our mask.
        Mask_t* mask = NULL;

        // If a mask file was given.
        if (eggsConfig -> maskFile != NULL) {

            // Read in whole VCF file from stdin.
            parse_vcf(replicate, inputStream);

            // Read in mask file.
            InputStream_t* maskInput = init_input_stream(eggsConfig -> maskFile);
            Replicate_t* maskReplicate = init_vcf_replicate(maskInput);
            parse_vcf(maskReplicate, maskInput);

            // Create mask.
            FourierCoefficients_t* fourierCoeff = init_fourier_coefficients(maskReplicate, eggsConfig -> numThreads);
            mask = create_fourier_mask(fourierCoeff, replicate -> numSamples, replicate -> numRecords, eggsConfig -> numThreads);

            // Destroy memory from mask replicate.
            destroy_fourier_coefficients(fourierCoeff);
            destroy_input_stream(maskInput);
            destroy_replicate(maskReplicate);

            // Apply fill to mask.
            if (eggsConfig -> fill > 0)
                apply_fill(replicate, mask, eggsConfig -> fill);

            print_replicate(replicate, mask, eggsConfig, fpOut);
        
        // If we want to create a random mask.
        } else if (eggsConfig -> randomMissing != NULL) {

            // Read in whole VCF file from stdin.
            parse_vcf(replicate, inputStream);
            // Create our random mask.
            mask = create_random_mask(replicate -> numSamples, replicate -> numRecords, eggsConfig -> meanMissing, eggsConfig -> stdMissing, eggsConfig -> numThreads);
            if (eggsConfig -> fill > 0)
                apply_fill(replicate, mask, eggsConfig -> fill);
            print_replicate(replicate, mask, eggsConfig, fpOut);

        // If no mask is created, then we can sequentially get a record and print it out.
        } else {
            Record_t* record = (Record_t*) calloc(1, sizeof(Record_t));
            record -> genotypes = (Genotype_t*) calloc(replicate -> numSamples, sizeof(Genotype_t));
            record -> numSamples = replicate -> numSamples;
            while (get_next_vcf_record(record, inputStream))
                print_record(record, NULL, eggsConfig, fpOut);
            destroy_record(record);
        }
        destroy_replicate(replicate);
        gzclose(fpOut);
        destroy_mask(mask);

    // Otherwise, parse as ms-style input.
    } else {
        Replicate_t* replicate = NULL;
        FourierCoefficients_t* fourierCoeff = NULL;

        // If a mask file was given, then get our Fourier coefficients.
        if (eggsConfig -> maskFile != NULL) {
            InputStream_t* maskInput = init_input_stream(eggsConfig -> maskFile);
            Replicate_t* maskReplicate = init_vcf_replicate(maskInput);
            parse_vcf(maskReplicate, maskInput);
            fourierCoeff = init_fourier_coefficients(maskReplicate, eggsConfig -> numThreads);
            destroy_input_stream(maskInput);
            destroy_replicate(maskReplicate);
        }

        // FOr each replicate in the ms-style input.
        while ((replicate = parse_ms(inputStream, eggsConfig -> length)) != NULL) {
            Mask_t* mask = NULL;
            // Generate Fourier mask.
            if (eggsConfig -> maskFile != NULL) {
                mask = create_fourier_mask(fourierCoeff, replicate -> numSamples, replicate -> numRecords, eggsConfig -> numThreads);
                if (eggsConfig -> fill > 0)
                    apply_fill(replicate, mask, eggsConfig -> fill);
            // Generate our random mask.
            } else if (eggsConfig -> randomMissing != NULL) {
                mask = create_random_mask(replicate -> numSamples, replicate -> numRecords, eggsConfig -> meanMissing, eggsConfig -> stdMissing, eggsConfig -> numThreads);
                if (eggsConfig -> fill > 0)
                    apply_fill(replicate, mask, eggsConfig -> fill);
            }
            gzFile fpOut;
            // If no basename given and only one replicate, then we print to stdout.
            if (ks_eof(inputStream -> fpIn) && eggsConfig -> outFile == NULL && numReps == 1) {
                fpOut = gzdopen(fileno(stdout), "w");
            // If output basename was given, then print replicates to basename.
            } else if (eggsConfig -> outFile != NULL) {
                kstring_t* outName = (kstring_t*) calloc(1, sizeof(kstring_t));
                ksprintf(outName, "%s_rep%d.vcf.gz", eggsConfig -> outFile, numReps);
                fpOut = gzopen(outName -> s, "w");
                free(outName -> s); free(outName);
            // If no basename was given use "rep".
            } else {
                kstring_t* outName = (kstring_t*) calloc(1, sizeof(kstring_t));
                ksprintf(outName, "rep%d.vcf.gz", numReps);
                fpOut = gzopen(outName -> s, "w");
                free(outName -> s); free(outName);
            }

            // Print out replicate.
            print_vcf_header(replicate, eggsConfig, fpOut);
            print_replicate(replicate, mask, eggsConfig, fpOut);

            gzclose(fpOut);
            destroy_replicate(replicate);
            destroy_mask(mask);
            numReps++;
        }
        destroy_fourier_coefficients(fourierCoeff);
    }

    destroy_eggs_configuration(eggsConfig);
    destroy_input_stream(inputStream);

    return 0;
}