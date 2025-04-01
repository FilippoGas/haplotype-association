# Filippo Gastaldello - 31/03/2025
#
# Build covariates file (.cov) used in Plink2 association analysis

library(tidyverse)
library(readxl)

# vcf input file
vcf <- snakemake@input[["vcf"]]
