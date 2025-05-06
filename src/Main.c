
#include "Interface.h"
#include "GenotypeParser.h"
#include "time.h"

typedef struct {

} Outputs

void print_vcf_header(Replicate_t* replicate, OutputStream_t* outputStream) {

}

int main(int argc, char* argv[]) {

    EggsConfig_t* eggsConfig = init_eggs_configuration(argc, argv);

    srand(time(NULL));

    InputStream_t* inputStream = init_input_stream(stdin);

    Replicate_t* replicate = init_vcf_replicate(inputStream);

    parse_vcf(replicate, inputStream);

    destroy_replicate(replicate);

    destroy_input_stream(inputStream);

    destroy_eggs_configuration(eggsConfig);

    return 0;
}