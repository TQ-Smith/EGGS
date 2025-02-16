# EGGS

Evolutionary Genotype Generalizer for Samples (EGGS) is a tool to remove phase, remove polarization for biallelic sites, and 
introduce missing data in genotypes. EGGS accepts variants in either VCF or ms-style replicates and produces
VCF output. For ms-style replicates, EGGS produces one VCF file per replicate.

## Building EGGS

The resulting executable will be in the **EGGS/bin** directory.

```
git clone https://github.com/TQ-Smith/EGGS.git 
cd EGGS
make
```

## Options

```
Usage: eggs [-upsc] [-l length] [-n number_of_samples] [-m prob_missing_allele] [-o out_basename]
Options:
   -l,--length      INT             Sets length of segment in number of base pairs for ms-replicates. 
                                        Default 1,000,000.
   -n,--numSamples  INT             Output the first n samples. 
                                        Default maximum possible per replicate.
   -m,--missing     DOUBLE          Genotypes are missing with supplied probability. 
                                        Default 0.
   -o,--out         STR             Basename for output files.
                                        Default "rep" prefix for multiple ms-style replicates.
   -u,--unphased                    Left and right genotypes are swapped with a probability of 50%.
   -p,--unpolarized                 Biallelic sites' ancestral and derived alleles with a probability of 50%.
   -s,--single                      Each sample contains one lineage. 
                                        ms-style input only. Ignores -u.
   -c,--compress                    Only used if -o is set. The resulting files are gzipped compressed.
```