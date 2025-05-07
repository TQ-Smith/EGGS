
#include "Interface.h"
#include "GenotypeParser.h"
#include "time.h"
#include "kstring.h"

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

void print_record(Record_t* record, EggsConfig_t* eggsConfig, gzFile fpOut) {
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

        if (eggsConfig -> hap) {
            record -> genotypes[i].right = MISSING;
            eggsConfig -> pseudohap = false;
        }

        if (eggsConfig -> unphase && !eggsConfig -> pseudohap) {
            record -> genotypes[i].isPhased = false;
            if (rand() < 0.5)
                SWAP(record -> genotypes[i].left, record -> genotypes[i].right, tempInt);
        }

        if (swapStates && record -> genotypes[i].left != MISSING)
            record -> genotypes[i].left ^= 1;
        if (swapStates && record -> genotypes[i].right != MISSING)
            record -> genotypes[i].right ^= 1;

        if (eggsConfig -> pseudohap && !eggsConfig -> hap) {
            record -> genotypes[i].isPhased = false;
            if (rand() < 0.5)
                record -> genotypes[i].right = record -> genotypes[i].left;
            else
                record -> genotypes[i].left = record -> genotypes[i].right;
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

void print_replicate(Replicate_t* replicate, EggsConfig_t* eggsConfig, gzFile fpOut) {
    Record_t* temp = replicate -> headRecord;
    for (int i = 0; i < replicate -> numRecords; i++) {
        print_record(temp, eggsConfig, fpOut);
        temp = temp -> nextRecord;
    }
}

int main(int argc, char* argv[]) {

    srand(time(NULL));

    EggsConfig_t* eggsConfig = init_eggs_configuration(argc, argv);
    if (eggsConfig == NULL)
        return -1;
    InputStream_t* inputStream = init_input_stream(stdin);
    int numReps = 1;

    ks_getuntil(inputStream -> fpIn, '\n', inputStream -> buffer, 0);
    if (strncmp(inputStream -> buffer -> s, "##fileformat=VCF", 16) == 0) {
        gzFile fpOut;
        if (eggsConfig -> outFile != NULL) {
            kstring_t* outName = calloc(1, sizeof(kstring_t));
            ksprintf(outName, "%s.vcf.gz", eggsConfig -> outFile);
            fpOut = gzopen(outName -> s, "w");
            free(outName -> s); free(outName);
        } else {
            fpOut = gzdopen(fileno(stdout), "w");
        }
        Replicate_t* replicate = init_vcf_replicate(inputStream);
        print_vcf_header(replicate, eggsConfig, fpOut);
        if (eggsConfig -> maskFile == NULL || eggsConfig -> randomMissing != NULL) {
            Record_t* record = calloc(1, sizeof(Record_t));
            record -> genotypes = calloc(replicate -> numSamples, sizeof(Genotype_t));
            record -> numSamples = replicate -> numSamples;
            while (get_next_vcf_record(record, inputStream))
                print_record(record, eggsConfig, fpOut);
            destroy_record(record);
        } else {
            parse_vcf(replicate, inputStream);
            print_replicate(replicate, eggsConfig, fpOut);
        }
        destroy_replicate(replicate);
        gzclose(fpOut);
    } else {
        Replicate_t* replicate = NULL;
        while ((replicate = parse_ms(inputStream, eggsConfig -> length)) != NULL) {
            gzFile fpOut;
            if (ks_eof(inputStream -> fpIn) && eggsConfig -> outFile == NULL && numReps == 1) {
                fpOut = gzdopen(fileno(stdout), "w");
            } else if (eggsConfig -> outFile != NULL) {
                kstring_t* outName = calloc(1, sizeof(kstring_t));
                ksprintf(outName, "%s_rep%d.vcf.gz", eggsConfig -> outFile, numReps);
                fpOut = gzopen(outName -> s, "w");
                free(outName -> s); free(outName);
            } else {
                kstring_t* outName = calloc(1, sizeof(kstring_t));
                ksprintf(outName, "rep%d.vcf.gz", numReps);
                fpOut = gzopen(outName -> s, "w");
                free(outName -> s); free(outName);
            }
            print_vcf_header(replicate, eggsConfig, fpOut);
            print_replicate(replicate, eggsConfig, fpOut);
            gzclose(fpOut);
            destroy_replicate(replicate);
            numReps++;
        }
    }

    destroy_eggs_configuration(eggsConfig);
    destroy_input_stream(inputStream);

    return 0;
}