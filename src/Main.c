
#include "Interface.h"
#include "GenotypeParser.h"
#include "time.h"

int main(int argc, char* argv[]) {

    EggsConfig_t* eggsConfig = init_eggs_configuration(argc, argv);

    srand(time(NULL));

    destroy_eggs_configuration(eggsConfig);

    return 0;
}