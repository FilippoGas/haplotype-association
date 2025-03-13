# Filippo Gastaldello - 12/02/2025
#
# Build covariates file (.cov) used in Plink2 association

library(tidyverse)
library(readxl)

# PATHS
vcf_file <- snakemake@input[["vcf"]]
pathway_alterations_file <- "/shares/CIBIO-Storage/BCG/scratch1/Resources/TCGA/dataTabs/somaticAlterations/10oncogenicPathways/1-s2.0-S0092867418303593-mmc4.xlsx"
tumor_subtype_file <- "/shares/CIBIO-Storage/BCG/scratch1/Resources/TCGA/dataTabs/Cancer_subtypes_list.tsv"
clinical_data_file <- "/shares/CIBIO-Storage/BCG/scratch1/Resources/TCGA/data/combined_study_clinical_data.tsv"
PCs_file <- "/shares/CIBIO-Storage/BCG/scratch1/Resources/TCGA/PCs/TCGA_Affy_all.eigenvec"
output_dir <- "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/data/TCGA_haplotype_association/BRCA_cancer_genes/"
tumor_type <- str_split_i(vcf_file, "_", 1)

# Phenotypes
tumor_subtypes <- read_delim(tumor_subtype_file, delim = "\t") %>%
                    filter(`TCGA PanCanAtlas Cancer Type Acronym` == tumor_type) %>% 
                    dplyr::select(COMMON, Subtype) %>% 
                    drop_na() %>%
                    pivot_wider(names_from = Subtype, values_from = Subtype)

# Widen subtype column in order to turn it into a binary variable since Plink2 
# accepts categorical variables, but with only 2 categories...
tumor_subtypes <- tumor_subtypes %>% column_to_rownames(var = "COMMON")
tumor_subtypes[!is.na(tumor_subtypes)] <- "1"
tumor_subtypes[is.na(tumor_subtypes)] <- "0"
tumor_subtypes <- tumor_subtypes %>% rownames_to_column(var = "COMMON")

pathway_alterations <- read_xlsx(path = pathway_alterations_file, sheet = "Pathway level")

# Covariates
age <- read_delim(clinical_data_file, delim = "\t") %>%
    filter(`TCGA PanCanAtlas Cancer Type Acronym` == tumor_type) %>% 
    dplyr::select(`Patient ID`, `Diagnosis Age`) %>% 
    dplyr::rename(patient = `Patient ID`, age = `Diagnosis Age`) %>% 
    unique()
sex <- read_delim(clinical_data_file, delim = "\t") %>%
    filter(`TCGA PanCanAtlas Cancer Type Acronym` == tumor_type) %>% 
    dplyr::select(`Patient ID`, `Sex`) %>% 
    dplyr::rename(patient = `Patient ID`) %>% 
    unique()
PCs <-  read_delim(PCs_file, delim = "\t") %>% 
        mutate(IID = str_sub(IID, 1, 16)) %>% 
    dplyr::select(-`#FID`)

# Retrieve sample list from newly generated vcf
samples_list <- sapply(colnames(as.data.frame(read_delim(file = vcf_file, delim = "\t")))[-(1:9)], function(x) gsub("\\.", "-", x))

# Build phenotype file and drop sample missing somatic information
pheno <- data.frame("sample" = as.character(samples_list)) %>% mutate(patient = str_sub(sample, 1, 12))
pheno <- pheno %>%  left_join(tumor_subtypes %>%
                                    dplyr::rename(patient = COMMON) %>%
                                    drop_na(),
                              by = "patient") %>%
                    left_join(pathway_alterations %>%
                                    mutate(patient = str_sub(SAMPLE_BARCODE, 1, 12)) %>% 
                                    dplyr::rename(Cell_Cycle = `Cell Cycle`,
                                                  RTK_RAS = `RTK RAS`),
                              by = "patient") %>% 
                    dplyr::select(-c(patient, SAMPLE_BARCODE)) %>% 
                    dplyr::rename(`#IID` = sample) %>% 
                    mutate(`#IID` = gsub("-",".", `#IID`))
pheno[is.na(pheno$Subtype), "Subtype"] <- "NONE"

# Build covariates file
cov <- data.frame("sample" = as.character(samples_list)) %>% mutate(patient = str_sub(sample, 1, 12), IID = str_sub(sample, 1, 16))
cov <- cov %>%  left_join(sex, by = "patient") %>% 
                left_join(age, by = "patient") %>% 
                left_join(PCs, by = "IID") %>% 
                dplyr::select(-c(patient, IID)) %>% 
                dplyr::rename(`#IID`= sample) %>% 
                mutate(`#IID` = gsub("-",".", `#IID`),
                       Sex = gsub("Female", "2", Sex),
                       Sex = gsub("Male", "1", Sex))


# Save .pheno
write_delim(pheno, file = snakemake@output[["pheno"]], delim = "\t", col_names = TRUE)
# Save .cov
write_delim(cov, file = snakemake@output[["cov"]], delim = "\t", col_names = TRUE)
