
#include "Missingness.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "time.h"
#include <pthread.h>

FourierCoefficients_t* init_fourier_coefficients(Replicate_t* replicate, int numThreads) {
    return NULL;
}

Mask_t* create_fourier_mask(FourierCoefficients_t* fourierCoeff, int numSamples, int numLoci, int numThreads) {
    return NULL;
}

Mask_t* create_random_mask(int numSamples, int numLoci, int numThreads) {
    return NULL;
}

void apply_fill(Replicate_t* replicate, Mask_t* mask, int fill) {

}

void destroy_mask(Mask_t* mask) {

}

void destroy_fourier_coefficients(int numSamples) {

}