# EGGS

Evolutionary Genotype Generalizer for Samples (EGGS) is a tool to remove phase, remove polarization for biallelic sites, and 
introduce missing alleles in genotypes. EGGS accepts variants in either VCF or ms-style replicates and produces
VCF output. For ms-style replicates, EGGS produces one VCF file per replicate. Output is gz-compressed.

## Building EGGS

The resulting executable will be in the **EGGS/bin** directory.

```
git clone https://github.com/TQ-Smith/EGGS.git 
cd EGGS
make
```

## Options

```
EGGS v1.0 February 2025
-----------------------

Written by T. Quinn Smith
Principal Investigator: Zachary A. Szpiech
The Pennsylvania State University

Usage: eggs [-upsch] [-l length] [-n number_of_samples] [-m prob_missing_allele] [-o out_basename]

Options:
   -l,--length      INT             Sets length of segment in number of base pairs for ms-replicates. 
                                        Default 1,000,000.
   -n,--numSamples  INT             Output the first n samples. 
                                        Default maximum possible per replicate.
   -m,--missing     DOUBLE          Genotypes are missing with supplied probability. 
                                        Default 0.
   -o,--out         STR             Basename for output files.
                                        Default "rep" prefix for multiple ms-style replicates.
   -u,--unphased                    Left and right genotypes are swapped with a probability of 0.5
   -p,--unpolarized                 Biallelic site alleles swapped with a probability of 0.5.
   -s,--single                      Each sample contains one lineage. 
                                        ms-style input only. Ignores -u.
   -h,--help                        Print help.
```