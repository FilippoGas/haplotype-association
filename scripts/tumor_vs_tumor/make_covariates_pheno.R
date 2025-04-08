# Filippo Gastaldello - 31/03/2025
#
# Build covariates file (.cov) used in Plink2 association analysis

library(tidyverse)
library(readxl)

####################
#   READ INPUTS    #
####################
vcf <- snakemake@input[["vcf"]]
clinical_data <- snakemake@input[["clinical_data"]]
PCs <- snakemake@input[["PCs"]]
tumor_list <- snakemake@params[["tumor_list"]]
# Get samples list from vcf
sample_list <- as.character(colnames(read_delim(vcf, delim = "\t", col_names = TRUE)))[-c(1:9)]

#############################
#   PREPARE PHENOTYPE FILE  #
#############################
pheno <- data.frame("sample" = sample_list) %>% mutate(`Patient ID` = str_sub(sample, 1, 12))
pheno <- pheno %>%  left_join(read_delim(clinical_data, delim = "\t", col_names = TRUE, col_types = "c") %>% 
                                  select(`Patient ID`, `TCGA PanCanAtlas Cancer Type Acronym`) %>%
                                  dplyr::rename(tumor = `TCGA PanCanAtlas Cancer Type Acronym`) %>% 
                                  filter(tumor %in% tumor_list) %>% 
                                  unique(),
                              by = "Patient ID") %>% 
                    select(-`Patient ID`) %>%
                    drop_na() %>% 
                    pivot_wider(values_from = "tumor", names_from = "tumor")
pheno <- pheno %>% column_to_rownames(var = "sample")
pheno[!is.na(pheno)] <- "1"
pheno[is.na(pheno)] <- "0"
pheno <- pheno %>% rownames_to_column(var = "#IID")
##############################
#   PREPARE COVARIATES FILE  #
##############################
cov <- data.frame("sample" = sample_list) %>% mutate(`Patient ID` = str_sub(sample, 1, 12))
cov <- cov %>% left_join(read_delim(clinical_data, delim = "\t", col_names = TRUE) %>% 
                                select(`Patient ID`, `Diagnosis Age`, Sex) %>% 
                                unique(),
                         by = "Patient ID") %>% 
                mutate(Sex = gsub("Male","1", Sex),
                       Sex = gsub("Female","2",Sex),
                       PCID = str_sub(sample, 1,19))
prcomp <- read_delim(PCs, delim = "\t", col_names = TRUE) %>% 
            select(-`#FID`) %>% 
            dplyr::rename("PCID" = "IID")
cov <- cov %>%  left_join(prcomp, by = "PCID") %>% 
                select(-c("Patient ID", "PCID")) %>% 
                dplyr::rename("#IID" = "sample",
                              "Age" = "Diagnosis Age")

write_tsv(cov, file = snakemake@output[["cov"]])
write_tsv(pheno, file = snakemake@output[["pheno"]])