
# File: plotMiss.R
# Date: 29 May 2025
# Author: T. Quinn Smith
# Principal Investigator: Zachary A. Szpiech
# Purpose: Visualize the missingness in a VCF file.

library(vcfR, quietly = TRUE, warn.conflicts = FALSE);
library(ggplot2, quietly = TRUE, warn.conflicts = FALSE);

# The name of VCF file is the only arg.
args <- commandArgs(trailingOnly = TRUE);
vcfFile <- args[2];
vcf <- read.vcfR(vcfFile, verbose = FALSE);

# Get the genotypes.
gt <- extract.gt(vcf, element = "GT");

# Test if a site is missing.
is_missing <- function(x) (is.na(x) | x == "./." | x == ".|.");

# Calculate missing proportion per locus
missing_prop <- apply(gt, 1, function(row) {mean(is_missing(row))});

# Create a dataframe with locus number and its associated proportion.
plot_data <- data.frame(Locus = 1:length(missing_prop), MissingProportion = missing_prop);

if (args[1] == "bar") {

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
    ggsave(paste(args[2], "_bar.png", sep = ""), plot = p, width = 15, height = 5, dpi = 300);

} else if (args[1] == "his") {
    breaks <- seq(0, 1, 0.005);
    p <- ggplot(plot_data, aes(x = MissingProportion)) +
        geom_histogram(aes(y = ..count.. / sum(..count..)), breaks=breaks, color = "blue") +
        ylab("Relative Frequency") +
        xlab("Proportion of Missing Samples Per Locus");

    ggsave(paste(args[2], "_his.png", sep = ""), plot = p, width = 10, height = 5, dpi = 300);
}