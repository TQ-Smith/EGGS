
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
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>

// Swap integer values.
#define SWAP(a, b, temp) { a = temp; a = b; b = temp; }

// Easy random uniform float with [0, 1)
#define rand() ((float) rand() / (float) (RAND_MAX))

// Convert eigen/ancestry map format to VCF file.
// Some code was taken from admixtools: https://github.com/DReichLab/AdmixTools/blob/master/src/geno.c
#define GENO_HEADER_SIZE 48
void convertEigenToVCF(char* genoFile, char* snpFile, char* indFile, gzFile fpOut) {
    gzFile ind = gzopen(indFile, "r");
    kstream_t* indStream = ks_init(ind);
    kstring_t* indBuffer = (kstring_t*) calloc(1, sizeof(kstring_t));

    // Print the VCF header information.
    gzprintf(fpOut, "##fileformat=VCFv4.2\n");
    gzprintf(fpOut, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    gzprintf(fpOut, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");   
    // Print the sample names in the header.
    int numSamplesInd = 0;
    while (true) {
        ks_getuntil(indStream, '\n', indBuffer, 0);
        if (ks_eof(indStream))
            break;
        char* indLine = strdup(indBuffer -> s);
        char* sampleName = strtok(indLine, " \t\n");
        gzprintf(fpOut, "\t%s", sampleName);
        free(indLine);
        numSamplesInd++;
    }
    gzprintf(fpOut, "\n");
    free(indBuffer -> s); free(indBuffer);
    ks_destroy(indStream);
    gzclose(ind);

    // For reading in geno file.
    bool transpose = false;
    char* genoMemStart = NULL;
    size_t genoMemLength = 0;
    unsigned char* genos = NULL;
    int numSnps = 0;
    int numSamples = 0;
    size_t rowSizeBytes;

    // Open geno file.
    struct stat sb;
    int fd = open(genoFile, O_RDONLY);
	if(fd < 0){
        fprintf(stderr, "Invalid genotype file. Exiting!\n");
		return;
	}

    // Map memory for file.
    fstat(fd, &sb);
	genoMemLength = (unsigned long) sb.st_size;
	genoMemStart = mmap(NULL, genoMemLength, PROT_READ, MAP_PRIVATE, fd, 0);
	if(genoMemStart == MAP_FAILED){
        close(fd);
        fprintf(stderr, "Could not read in entire genotype file. Exiting!\n");
		return;
	}
	close(fd);

    char header[GENO_HEADER_SIZE+1];
    header[GENO_HEADER_SIZE] = 0;
    char matrixType[GENO_HEADER_SIZE];
    strncpy(header, genoMemStart, GENO_HEADER_SIZE);
    sscanf(header, "%s %d %d", matrixType, &numSamples, &numSnps);
    // Check number of samples.
    if (numSamples != numSamplesInd) {
        munmap(genoMemStart, genoMemLength);
        fprintf(stderr, "Unequal number of samples between geno and ind files. Exiting!\n");
        return;
    }

    size_t headerSize = GENO_HEADER_SIZE;
    if(strncmp(matrixType, "GENO", 4) == 0){
		transpose = false;
        size_t neededBytes = (size_t) ceil(numSamples / 4.0);
        rowSizeBytes = (neededBytes > GENO_HEADER_SIZE) ? neededBytes : GENO_HEADER_SIZE;
        headerSize = rowSizeBytes;
	} else if (strncmp(matrixType, "TGENO", 5) == 0){
		transpose = true;
		rowSizeBytes = (size_t) ceil(numSnps / 4.0);
	} else{
        munmap(genoMemStart, genoMemLength);
		fprintf(stderr, "Invalid file type. Exiting!\n");
		return;
	}
    genos = (unsigned char*) genoMemStart + headerSize;


    // Now, we parse each site.
    gzFile snp = gzopen(snpFile, "r");
    kstream_t* snpStream = ks_init(snp);
    kstring_t* snpBuffer = (kstring_t*) calloc(1, sizeof(kstring_t));
    int numSnpsSnps = 0;
    while (true) {
        ks_getuntil(snpStream, '\n', snpBuffer, 0);
        if (ks_eof(snpStream) || numSnpsSnps > numSnps)
            break;
        // Get all info of the SNP.
        char* snpLine = strdup(snpBuffer -> s);
        char* token = strtok(snpLine, " \t\n");
        char* id = strdup(token); 
        token = strtok(NULL, " \t\n");
        char* chrom = strdup(token);
        token = strtok(NULL, " \t\n");
        token = strtok(NULL, " \t\n");
        char* position = strdup(token);
        token = strtok(NULL, " \t\n");
        char* ref = strdup(token);
        token = strtok(NULL, " \t\n");
        char* alt = strdup(token);

        // Print the site information.
        gzprintf(fpOut, "%s\t%s\t%s\t%s\t%s\t.\t.\t.\tGT", chrom, position, id, ref, alt);
        free(snpLine); free(id); free(chrom); free(position); free(ref); free(alt);

        // Print each sample's genotype.
        size_t byteIndex; 
        int bitShift;
        for (int i = 0; i < numSamples; i++) {
            if (!transpose) {
                byteIndex = i / 4; 
                bitShift = 6 - (i % 4) * 2;
                byteIndex += (rowSizeBytes * numSnpsSnps);
            } else {
                byteIndex = numSnpsSnps / 4; 
                bitShift = 6 - (numSnpsSnps % 4) * 2;
                byteIndex += (rowSizeBytes * i);
            }
            unsigned char genoByte = genos[byteIndex];
            switch (((genoByte >> bitShift) & 0x3)) {
                case 0: gzprintf(fpOut, "\t0/0"); break;
                case 1: gzprintf(fpOut, "\t0/1"); break;
                case 2: gzprintf(fpOut, "\t1/1"); break;
                case 3: gzprintf(fpOut, "\t./."); break;
            }

        }
        gzprintf(fpOut, "\n");
        numSnpsSnps++;
    }

    if (numSnpsSnps != numSnps)
        fprintf(stderr, "Unequal number of snps between geno and snp files. Exiting!\n");

    munmap(genoMemStart, genoMemLength);
    free(snpBuffer -> s); free(snpBuffer);
    ks_destroy(snpStream); 
    gzclose(snp);
}

void summary(Replicate_t* replicate, InputStream_t* inputStream, char* outName) {
    // Standard error by default.
    FILE* indOut = stdout;
    FILE* lociOut = stdout;
    if (outName != NULL) {
        kstring_t* out = calloc(1, sizeof(kstring_t));
        ksprintf(out, "%s.ind.tsv", outName);
        indOut = fopen(out -> s, "w");
        free(out -> s); free(out);
        kstring_t* out2 = calloc(1, sizeof(kstring_t));
        ksprintf(out2, "%s.loci.tsv", outName);
        lociOut = fopen(out2 -> s, "w");
        free(out2 -> s); free(out2);
    }

    int* sampleProportions = calloc(replicate -> numSamples, sizeof(int));

    Record_t* record = (Record_t*) calloc(1, sizeof(Record_t));
    record -> genotypes = (Genotype_t*) calloc(replicate -> numSamples, sizeof(Genotype_t));
    record -> numSamples = replicate -> numSamples;
    
    // Print loci file.
    fprintf(lociOut, "CHROM\tPOS\tPROP_MISSING\n");
    int numRecords = 0;
    while (get_next_vcf_record(record, inputStream, false)) {
        int numMissing = 0;
        for (int i = 0; i < replicate -> numSamples; i++) {
            if (record -> genotypes[i].left == MISSING && record -> genotypes[i].right == MISSING) {
                sampleProportions[i] += 1;
                numMissing++;
            }
        }
        fprintf(lociOut, "%s\t%d\t%lf\n", record -> chrom, record -> position, numMissing / (double) replicate -> numSamples);
        numRecords++;
    }

    if (outName == NULL) fprintf(stdout, "\n");

    // Print sample files.
    fprintf(indOut, "SAMPLE\tPROP_MISSING\n");
    for (int i = 0; i < replicate -> numSamples; i++)
        fprintf(indOut, "%s\t%lf\n", replicate -> sampleNames[i], sampleProportions[i] / (double) numRecords);

    destroy_record(record);
    free(sampleProportions);
    fclose(indOut);
    fclose(lociOut);
}

// Print a simple VCF header from a replicate.
void print_vcf_header(Replicate_t* replicate, EggsConfig_t* eggsConfig, gzFile fpOut) {
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
//  int* mask -> If index contains MISSING, then the sample's genotype is masked out.
//  EggsConfig_t* eggsConfig -> Defines the logic of the record manipulation.
//  gzFile fpOut -> The output stream.
// Returns: void.
void print_record(Record_t* record, int* mask, EggsConfig_t* eggsConfig, gzFile fpOut) {
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
    if (eggsConfig -> keep)
        gzprintf(fpOut, "\t.\t.\t%s\tGT", record -> info);
    else
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
        if (mask != NULL && mask[i] == MISSING) {
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

            // If sequencing error.
            if (!eggsConfig -> pseudohap && eggsConfig -> probTransition == 0 && record -> numAlleles == 2 && eggsConfig -> seqerr > 0) {
                if (record -> genotypes[i].left != MISSING && rand() < eggsConfig -> seqerr)
                    record -> genotypes[i].left ^= 1;
                if (record -> genotypes[i].right != MISSING && rand() < eggsConfig -> seqerr)
                    record -> genotypes[i].right ^= 1;
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

            // If missing, use ancestral. If multiallelic use alternative.
            if (temp -> genotypes[i].left > 1 || temp -> genotypes[i].left == MISSING)
                temp -> genotypes[i].left = 1;
            if (temp -> genotypes[i].right > 1 || temp -> genotypes[i].right == MISSING)
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
        for (Record_t* temp = replicate -> headRecord; temp != NULL; temp = temp -> nextRecord)
            gzprintf(fpOut, "%d", temp -> genotypes[i].right);
        gzprintf(fpOut, "\n");
    }
}

// Print the list of records in a replicate.
// Accepts:
//  Replicate_t* replicate -> The records to print.
//  MissingDistribution_t* mask -> Used to create the mask.
//  EggsConfig_t* eggsConfig -> The user parameters.
//  gzFile fpOut -> The output stream.
// Returns: void.
void print_replicate(Replicate_t* replicate, MissingDistribution_t* mask, EggsConfig_t* eggsConfig, gzFile fpOut) {
    if (!eggsConfig -> msOutput) {
        int* site = NULL;
        if (mask != NULL)
            site = calloc(replicate -> numSamples, sizeof(int));
        Record_t* temp = replicate -> headRecord;
        for (int i = 0; i < replicate -> numRecords; i++) {
            // If we are choosing randomly from the realized distribution.
            if (eggsConfig -> randomMissing != NULL)
                create_random_mask(mask, replicate -> numSamples, 0, replicate -> numRecords, eggsConfig -> meanMissing, eggsConfig -> stdMissing, site);
            // If we are using beta or EGGS methods.
            if (eggsConfig -> betaMissing != NULL || eggsConfig -> maskFile != NULL)
                create_random_mask(mask, replicate -> numSamples, i, replicate -> numRecords, eggsConfig -> meanMissing, eggsConfig -> stdMissing, site);
            print_record(temp, site, eggsConfig, fpOut);
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
    Replicate_t* replicate = init_vcf_replicate(inputStream, true);
    
    MissingDistribution_t* dis = init_missing_distribution(replicate, inputStream);
    
    // Calculate and set mean.
    double mean = 0;
    for (int i = 0; i < dis -> numRecords; i++)
        mean += dis -> proportions[i];
    *mu = mean / (double) dis -> numRecords;

    // Calculate and set variance.
    double var = 0;
    for (int i = 0; i < dis -> numRecords; i++)
        var += (dis -> proportions[i] - *mu) * (dis -> proportions[i] - *mu);
    *sigma = sqrt(var / (dis -> numRecords - 1.0));

    destroy_missing_distribution(dis);
    destroy_replicate(replicate);
}

int main(int argc, char* argv[]) {

    srand(time(NULL));

    // Get CLI configuration.
    EggsConfig_t* eggsConfig = init_eggs_configuration(argc, argv);
    
    if (eggsConfig == NULL)
        return -1;
    
    // If -e was set, convert file and exit.
    if (eggsConfig -> eigenFiles != NULL) {
        char* fileNames = strdup(eggsConfig -> eigenFiles);
        char* genoFile = strtok(fileNames, ",");
        char* snpFile = strtok(NULL, ",");
        char* indFile = strtok(NULL, ",");
        // Open output stream.
        gzFile fpOut;
        if (eggsConfig -> outFile != NULL) {
            kstring_t* outName = (kstring_t*) calloc(1, sizeof(kstring_t));
            ksprintf(outName, "%s.vcf.gz", eggsConfig -> outFile);
            fpOut = gzopen(outName -> s, "w");
            free(outName -> s); free(outName);
        } else {
            fpOut = gzdopen(fileno(stdout), "w");
        }
        // Convert and print.
        convertEigenToVCF(genoFile, snpFile, indFile, fpOut);
        // Exit.
        free(fileNames);
        destroy_eggs_configuration(eggsConfig);
        gzclose(fpOut);
        return 1;
    }
    
    // If a VCF file was given for random missingness, calculate mean and standard deviation per site.
    if (eggsConfig -> betaMissing != NULL && eggsConfig -> meanMissing == -1 && eggsConfig -> stdMissing == -1) {
        InputStream_t* betaStream = init_input_stream(eggsConfig -> betaMissing);
        get_mu_sigma(betaStream, &(eggsConfig -> meanMissing), &(eggsConfig -> stdMissing));
        destroy_input_stream(betaStream);
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
    
    // Read in distribution if supplied.
    MissingDistribution_t* mask = NULL;
    InputStream_t* maskInput = NULL;
    Replicate_t* missingReplicate = NULL;
    if (eggsConfig -> maskFile != NULL) {
        maskInput = init_input_stream(eggsConfig -> maskFile);
        missingReplicate = init_vcf_replicate(maskInput, true);
        mask = init_missing_distribution(missingReplicate, maskInput);
        destroy_input_stream(maskInput);
        destroy_replicate(missingReplicate);
    } else if (eggsConfig -> randomMissing != NULL) {
        maskInput = init_input_stream(eggsConfig -> randomMissing);
        missingReplicate = init_vcf_replicate(maskInput, true);
        mask = init_missing_distribution(missingReplicate, maskInput);
        destroy_input_stream(maskInput);
        destroy_replicate(missingReplicate);
    } else if (eggsConfig -> betaMissing != NULL) {
        mask = calloc(1, sizeof(MissingDistribution_t));
        // Beta distribution just needs rng.
        gsl_rng_env_setup();
        const gsl_rng_type* T = gsl_rng_default;
        gsl_rng* r = gsl_rng_alloc(T);
        gsl_rng_set(r, time(NULL));
        mask -> r = r;
    }
    
    // If VCF.
    if (strncmp(inputStream -> buffer -> s, "##fileformat=VCF", 16) == 0) {

        // If verbose summary statistics only.
        if (eggsConfig -> verbose) {
            Replicate_t* replicate = init_vcf_replicate(inputStream, true);
            summary(replicate, inputStream, eggsConfig -> outFile);
            destroy_replicate(replicate);
            destroy_input_stream(inputStream);
            destroy_eggs_configuration(eggsConfig);
            return 0;
        }

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

        // Eat meta info lines until line starting with #CHROM.
        gzprintf(fpOut, "##fileformat=VCFv4.2\n");
        int dret;
        do {
            ks_getuntil(inputStream -> fpIn, '\n', inputStream -> buffer, &dret);
            // Echo the contig info.
            if (strncmp(inputStream -> buffer -> s, "##CONTIG", 8) == 0)
                gzprintf(fpOut, "%s\n", inputStream -> buffer -> s);
            // Echo info tag if flag set.
            if (eggsConfig -> keep && strncmp(inputStream -> buffer -> s, "##INFO", 6) == 0)
                gzprintf(fpOut, "%s\n", inputStream -> buffer -> s);
        } while (strncmp(inputStream -> buffer -> s, "#C", 2) != 0);

        Replicate_t* replicate = init_vcf_replicate(inputStream, false);
        parse_vcf(replicate, inputStream, eggsConfig -> keep);
        print_vcf_header(replicate, eggsConfig, fpOut);
        print_replicate(replicate, mask, eggsConfig, fpOut);
        destroy_replicate(replicate);
        gzclose(fpOut);

    // Otherwise, parse as ms-style input.
    } else {
        
        Replicate_t* replicate = NULL;

        // For each replicate in the ms-style input.
        int length;
        if (eggsConfig -> length <= 0 )
            length = 1000000;
        else 
            length = eggsConfig -> length;
        
        while ((replicate = parse_ms(inputStream, length, eggsConfig -> hap)) != NULL) {
            
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
            if (!eggsConfig -> msOutput)
                print_vcf_header(replicate, eggsConfig, fpOut);
            print_replicate(replicate, mask, eggsConfig, fpOut);

            gzclose(fpOut);
            destroy_replicate(replicate);
            numReps++;
        }
    }

    destroy_eggs_configuration(eggsConfig);
    destroy_input_stream(inputStream);
    destroy_missing_distribution(mask);

    return 0;
}