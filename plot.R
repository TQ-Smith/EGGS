library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

df <- read.delim(args[1], header = TRUE, sep = "\t")

p <- ggplot(df, aes(x = RECORD, y = PROPORTION)) +
  geom_line(color = "blue", size = 1.2) +
  labs(title = "Proportion of Samples with Missing Genotypes at Each Locus",
       x = "Locus",
       y = "Proportion of Samples with Missing Genotypes") +
  theme_classic()
ggsave(paste(args[1], ".png", sep=""), plot = p, width = 10, height = 4, dpi = 300)
