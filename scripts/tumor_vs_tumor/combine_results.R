# Filippo Gastaldello - 10/04/2025
#
# Combine Plink results of different chromosomes for all tumor/model couples

library(tidyverse)

# get Plink results for chr1 of a given model/phenotype 
chr1_result <- snakemake@input[["plink_results"]]
# make list of all chromosome's inputs
chr1_result <- str_split(chr1_result, pattern = "chr1")
plink_results_list <- c()
for (chr in seq(1,22)) {
    plink_results_list <- c(plink_results_list, as.character(paste0(chr1_result[[1]][1],"chr",chr,chr1_result[[1]][2])))
}
# Open and concat results from all chromosomes
plink_results <- data_frame()
for (file in plink_results_list) {
    plink_results <- bind_rows(plink_results, read_tsv(file))
}
write_tsv(plink_results, file = snakemake@output[[1]])