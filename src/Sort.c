
// File: Sort.c
// Date: 20 May 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Sorting methods.

void swapDouble(double* a, double* b) {
    double temp = *a;
    *a = *b;
    *b = temp;
}

void swapInt(int* a, int* b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

int partition_bins(double* avgPower, int* sortedBins, int low, int high) {
    double pivot = avgPower[high]; 
    int i = (low - 1);
    for (int j = low; j < high; j++) {
        if (avgPower[j] >= pivot) {
            i++;
            swapDouble(&avgPower[i], &avgPower[j]);
            swapInt(&sortedBins[i], &sortedBins[j]);
        }
    }
    swapDouble(&avgPower[i + 1], &avgPower[high]);
    swapInt(&sortedBins[i + 1], &sortedBins[high]);
    return (i + 1);
}

void quick_sort_bins(double* avgPower, int* sortedBins, int low, int high) {
    if (low < high) {
        int pi = partition_bins(avgPower, sortedBins, low, high);
        quick_sort_bins(avgPower, sortedBins, low, pi - 1);
        quick_sort_bins(avgPower, sortedBins, pi + 1, high);
    }
}

int partition_indices(int* binIndices, int low, int high) {
    int pivot = binIndices[high]; 
    int i = (low - 1);
    for (int j = low; j < high; j++) {
        if (binIndices[j] <= pivot) {
            i++;
            swapInt(&binIndices[i], &binIndices[j]);
        }
    }
    swapInt(&binIndices[i + 1], &binIndices[high]);
    return (i + 1);
}

void quick_sort_indices(int* binIndices, int low, int high) {
    if (low < high) {
        int pi = partition_indices(binIndices, low, high);
        quick_sort_indices(binIndices, low, pi - 1);
        quick_sort_indices(binIndices, pi + 1, high);
    }
}