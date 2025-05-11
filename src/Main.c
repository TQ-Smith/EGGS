
#include "Interface.h"
#include "GenotypeParser.h"
#include "time.h"
#include "kstring.h"
#include "Missingness.h"

#define SWAP(a, b, temp) { a = temp; a = b; b = temp; }

#define rand() ((float) rand() / (float) (RAND_MAX))

void print_vcf_header(Replicate_t* replicate, EggsConfig_t* eggsConfig, gzFile fpOut) {
    gzprintf(fpOut, "##fileformat=VCFv4.2\n"); 
    gzprintf(fpOut, "##eggsCommand=%s\n", eggsConfig -> command);
    gzprintf(fpOut, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"); 
    if (replicate -> sampleNames == NULL)
        for (int i = 0; i < replicate -> numSamples; i++)
            gzprintf(fpOut, "\ts%d", i); 
    else
        for (int i = 0; i < replicate -> numSamples; i++)
            gzprintf(fpOut, "\t%s", replicate -> sampleNames[i]);
    gzprintf(fpOut, "\n");
}

void print_record(Record_t* record, int* mask, EggsConfig_t* eggsConfig, gzFile fpOut) {
    if (record -> chrom == NULL) 
        gzprintf(fpOut, "chr1");
    else 
        gzprintf(fpOut, "%s", record -> chrom);
    gzprintf(fpOut, "\t%d\t.", record -> position);
    if (record -> ref == NULL)
        gzprintf(fpOut, "\tA");
    else 
        gzprintf(fpOut, "\t%s", record -> ref);
    if (record -> alts == NULL)
        gzprintf(fpOut, "\tT");
    else 
        gzprintf(fpOut, "\t%s", record -> alts);
    gzprintf(fpOut, "\t.\t.\t.\t.");

    int tempInt = 0;
    bool swapStates = false;
    if (eggsConfig -> unpolarize && record -> numAlleles == 2 && rand() < 0.5)
        swapStates = true;
    for (int i = 0; i < record -> numSamples; i++) {

        if (mask != NULL && mask[i] == MISSING) {
            record -> genotypes[i].isPhased = false;
            record -> genotypes[i].left = MISSING;
            record -> genotypes[i].right = MISSING;
        } else {

            if (eggsConfig -> unphase && !eggsConfig -> pseudohap) {
                record -> genotypes[i].isPhased = false;
                if (rand() < 0.5)
                    SWAP(record -> genotypes[i].left, record -> genotypes[i].right, tempInt);
            }

            if (swapStates && record -> genotypes[i].left != MISSING)
                record -> genotypes[i].left ^= 1;
            if (swapStates && record -> genotypes[i].right != MISSING)
                record -> genotypes[i].right ^= 1;

            if (eggsConfig -> pseudohap) {
                record -> genotypes[i].isPhased = false;
                if (rand() < 0.5)
                    record -> genotypes[i].right = record -> genotypes[i].left;
                else
                    record -> genotypes[i].left = record -> genotypes[i].right;
            }

        }

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

void print_replicate(Replicate_t* replicate, Mask_t* mask, EggsConfig_t* eggsConfig, gzFile fpOut) {
    Record_t* temp = replicate -> headRecord;
    for (int i = 0; i < replicate -> numRecords; i++) {
        if (mask == NULL)
            print_record(temp, NULL, eggsConfig, fpOut);
        else 
            print_record(temp, mask -> missing[i], eggsConfig, fpOut);
        temp = temp -> nextRecord;
    }
}

int main(int argc, char* argv[]) {

    srand(time(NULL));

    EggsConfig_t* eggsConfig = init_eggs_configuration(argc, argv);
    if (eggsConfig == NULL)
        return -1;
    InputStream_t* inputStream = init_input_stream(NULL);
    int numReps = 1;

    ks_getuntil(inputStream -> fpIn, '\n', inputStream -> buffer, 0);
    if (strncmp(inputStream -> buffer -> s, "##fileformat=VCF", 16) == 0) {
        gzFile fpOut;
        if (eggsConfig -> outFile != NULL) {
            kstring_t* outName = (kstring_t*) calloc(1, sizeof(kstring_t));
            ksprintf(outName, "%s.vcf.gz", eggsConfig -> outFile);
            fpOut = gzopen(outName -> s, "w");
            free(outName -> s); free(outName);
        } else {
            fpOut = gzdopen(fileno(stdout), "w");
        }
        Replicate_t* replicate = init_vcf_replicate(inputStream);
        Mask_t* mask = NULL;
        print_vcf_header(replicate, eggsConfig, fpOut);
        if (eggsConfig -> maskFile == NULL && eggsConfig -> randomMissing == NULL) {
            Record_t* record = (Record_t*) calloc(1, sizeof(Record_t));
            record -> genotypes = (Genotype_t*) calloc(replicate -> numSamples, sizeof(Genotype_t));
            record -> numSamples = replicate -> numSamples;
            while (get_next_vcf_record(record, inputStream))
                print_record(record, NULL, eggsConfig, fpOut);
            destroy_record(record);
        } else {
            parse_vcf(replicate, inputStream);
            if (eggsConfig -> maskFile != NULL) {
                InputStream_t* maskInput = init_input_stream(eggsConfig -> maskFile);
                Replicate_t* maskReplicate = init_vcf_replicate(maskInput);
                parse_vcf(maskReplicate, maskInput);
                FourierCoefficients_t* fourierCoeff = init_fourier_coefficients(maskReplicate);
                mask = create_fourier_mask(fourierCoeff, replicate -> numSamples, replicate -> numRecords);
                for (int i = 0; i < replicate -> numRecords; i++) {
                    for (int j = 0; j < replicate -> numSamples; j++) {
                        gzprintf(fpOut, "%d\t", mask ->missing[j][i]);
                    }
                    gzprintf(fpOut, "\n");
                }
                destroy_fourier_coefficients(fourierCoeff);
                destroy_input_stream(maskInput);
                destroy_replicate(maskReplicate);
                if (eggsConfig -> fill > 0)
                    apply_fill(replicate, mask, eggsConfig -> fill);
            } else if (eggsConfig -> randomMissing != NULL) {
                mask = create_random_mask(replicate -> numSamples, replicate -> numRecords, eggsConfig -> meanMissing, eggsConfig -> stdMissing);
                if (eggsConfig -> fill > 0)
                    apply_fill(replicate, mask, eggsConfig -> fill);
            }
            print_replicate(replicate, mask, eggsConfig, fpOut);
        }
        destroy_replicate(replicate);
        gzclose(fpOut);
        destroy_mask(mask);
    } else {
        Replicate_t* replicate = NULL;
        FourierCoefficients_t* fourierCoeff = NULL;
        if (eggsConfig -> maskFile != NULL) {
            InputStream_t* maskInput = init_input_stream(eggsConfig -> maskFile);
            Replicate_t* maskReplicate = init_vcf_replicate(maskInput);
            parse_vcf(maskReplicate, maskInput);
            fourierCoeff = init_fourier_coefficients(maskReplicate);
            destroy_input_stream(maskInput);
            destroy_replicate(maskReplicate);
        }
        while ((replicate = parse_ms(inputStream, eggsConfig -> length)) != NULL) {
            Mask_t* mask = NULL;
            if (eggsConfig -> maskFile != NULL) {
                mask = create_fourier_mask(fourierCoeff, replicate -> numSamples, replicate -> numRecords);
                if (eggsConfig -> fill > 0)
                    apply_fill(replicate, mask, eggsConfig -> fill);
            } else if (eggsConfig -> randomMissing != NULL) {
                mask = create_random_mask(replicate -> numSamples, replicate -> numRecords, eggsConfig -> meanMissing, eggsConfig -> stdMissing);
                if (eggsConfig -> fill > 0)
                    apply_fill(replicate, mask, eggsConfig -> fill);
            }
            gzFile fpOut;
            if (ks_eof(inputStream -> fpIn) && eggsConfig -> outFile == NULL && numReps == 1) {
                fpOut = gzdopen(fileno(stdout), "w");
            } else if (eggsConfig -> outFile != NULL) {
                kstring_t* outName = (kstring_t*) calloc(1, sizeof(kstring_t));
                ksprintf(outName, "%s_rep%d.vcf.gz", eggsConfig -> outFile, numReps);
                fpOut = gzopen(outName -> s, "w");
                free(outName -> s); free(outName);
            } else {
                kstring_t* outName = (kstring_t*) calloc(1, sizeof(kstring_t));
                ksprintf(outName, "rep%d.vcf.gz", numReps);
                fpOut = gzopen(outName -> s, "w");
                free(outName -> s); free(outName);
            }
            print_vcf_header(replicate, eggsConfig, fpOut);
            print_replicate(replicate, mask, eggsConfig, fpOut);
            gzclose(fpOut);
            destroy_replicate(replicate);
            destroy_fourier_coefficients(fourierCoeff);
            destroy_mask(mask);
            numReps++;
        }
    }

    destroy_eggs_configuration(eggsConfig);
    destroy_input_stream(inputStream);

    return 0;
}