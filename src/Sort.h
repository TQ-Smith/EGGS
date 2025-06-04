
// File: Sort.h
// Date: 20 May 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Sorting methods.

// These methods are repetitive but straightforward.

// Sort avgPower from greatest to least. Move sortedBins around along with avgPower.
void quick_sort_bins(double* avgPower, int* sortedBins, int low, int high);

// Traditional quick sort of integers from least to greatest.
void quick_sort_indices(int* binIndices, int low, int high);