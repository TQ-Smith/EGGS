
# File: plotMiss.R
# Date: 29 May 2025
# Author: T. Quinn Smith
# Principal Investigator: Zachary A. Szpiech
# Purpose: Plot the proportion of missing samples (./. or .|.) at each VCF record.

library(vcfR, quietly = TRUE, warn.conflicts = FALSE);
library(ggplot2, quietly = TRUE, warn.conflicts = FALSE);

# The name of VCF file is the only arg.
args <- commandArgs(trailingOnly = TRUE);
vcfFile <- args[1];
vcf <- read.vcfR(vcfFile, verbose = FALSE);

# Get the genotypes.
gt <- extract.gt(vcf, element = "GT");

# Test if a site is missing.
is_missing <- function(x) (is.na(x) | x == "./." | x == ".|.");

# Calculate missing proportion per locus
missing_prop <- apply(gt, 1, function(row) {mean(is_missing(row))});

# Create a dataframe with locus number and its associated proportion.
plot_data <- data.frame(Locus = 1:length(missing_prop), MissingProportion = missing_prop);

# Create plot from dataframe.
p <- ggplot(plot_data, aes(x = factor(Locus), y = MissingProportion)) +
    geom_bar(stat = "identity", color = "blue") +
    labs(
        title = "Proportion of Missing Samples Per Locus",
        x = "Locus",
        y = "Proportion"
    ) +
    theme( 
        text = element_text(size = 12),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    ) +
    theme(axis.text.x = element_blank());

# Save our plot.
ggsave(paste(args[1], "_missing.png", sep = ""), plot = p, width = 15, height = 5, dpi = 300);