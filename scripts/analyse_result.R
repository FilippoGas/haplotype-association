# Filippo Gastaldello - 19/02/2025
#
# Analyse result from GWAS

library(qqman)

gwas_result <- read_delim(file = "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/data/TCGA_haplotype_association/BRCA_cancer_genes/results/BRCA_cancer_genes_ADDITIVE.BRCA_Her2.glm.logistic.hybrid",
                          delim = "\t")
gwas_result <- gwas_result %>% filter(TEST == "ADD") %>% drop_na()
manhattan(gwas_result, chr = "#CHROM" ,bp = "POS", snp = "ID")


