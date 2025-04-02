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
sample_list <- colnames(as_data_frame(read_delim(vcf, delim = "\t", col_names = TRUE))[-c(1:9)])

#############################
#   PREPARE PHENOTYPE FILE  #
#############################
tumor_type <- read_delim(clinical_data, delim = "\t", col_names = TRUE) %>% 
            select(`Patient ID`, `TCGA PanCanAtlas Cancer Type Acronym`) %>%
            dplyr::rename(tumor = `TCGA PanCanAtlas Cancer Type Acronym`) %>% 
            filter(tumor %in% tumor_list) %>% 
            pivot_wider(names_from = tumor, values_from = tumor)
tumor_type <- tumor_type %>% column_to_rownames(var = "Patient ID")
tumor_type[!is.na(tumor_type)] <- "1" 
tumor_type[is.na(tumor_type)] <- "0" 
tumor_type <- tumor_type %>% rownames_to_column(var = "Patient ID")
pheno <- data.frame("sample" = sample_list) %>% mutate(`Patient ID` = str_sub(sample, 1, 12))
pheno <- pheno %>%  left_join(tumor_type, by = "Patient ID") %>% 
                    select(-`Patient ID`) %>% 
                    dplyr::rename(`#IID` = "sample")

##############################
#   PREPARE COVARIATES FILE  #
##############################
cov <- data.frame("sample" = sample_list) %>% mutate(`Patient ID` = str_sub(sample, 1, 12))
cov <- cov %>% left_join(read_delim(clinical_data, delim = "\t", col_names = TRUE) %>% 
                            select(`Patient ID`, `Diagnosis Age`, Sex),
                         by = "Patient ID") %>% 
                mutate(Sex = gsub("Male","1", Sex),
                       Sex = gsub("Female","2",Sex),
                       PCID = str_sub(sample, 1,19))
prcomp <- read_delim(PCs, delim = "\t", col_names = TRUE) %>% 
            select(-`#FID`) %>% 
            dplyr::rename("PCID" = "IID")
cov <- cov %>%  left_join(prcomp, by = "PCID") %>% 
                select(-c("Patient ID", "PCID"))

write_tsv(cov, file = snakemake@output[["cov"]])
write_tsv(pheno, file = snakemake@output[["pheno"]])