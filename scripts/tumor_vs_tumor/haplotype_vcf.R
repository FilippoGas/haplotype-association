# Filippo Gastaldello - 13/03/2025
#
# In order to conduct tumor vs tumor association in Plink, we need a vcf file 
# containing all samples from all tumors together, indicating their genotype for 
# each protein coding transcript.

library(tidyverse)
library(parallel)

# Get location of file
snakefile_abs_path <- gsub("[\r\n]", "", snakemake@params[["abs_path"]])
# Maximum number of ALT alleles (Plink2 constrain)
ALT_limit = 254
# Get number of cores to use from params
cores <- snakemake@params[["cores"]]
# Get chromosome name from expected output name
chr <- str_split_i(str_split_i(basename(snakemake@output[["vcf"]]), "_", 1), "chr", 2)

################################
# READ INPUTS FROM CONFIG.YAML #
################################
config <- yaml::read_yaml(paste0(snakefile_abs_path,'/config.yaml'))
# List of tumor types to take into consideration
tumor_list <- config$cancer_types
# wt sequence of proteins of interest, only keep those with a valid sequence
load(paste0(config$wt_transcripts, "/wt_aa.RData"))
sequences_aa <- column_to_rownames(sequences_aa, var = "ensembl_transcript_id")
sequences_aa <- sequences_aa %>% filter(sequence != "Sequence unavailable")
# transcript annotations (actually only need transcript start and chromosome to
# order them in the vcf file)
load(paste0(config$wt_transcripts, "/protein_coding_transcripts.RData"))
transcript_positions <- protein_coding_transcripts %>%
                            filter(ensembl_transcript_id %in% rownames(sequences_aa),
                                   chromosome_name == chr) %>%
                            select(ensembl_transcript_id, transcript_start, chromosome_name) %>% 
                            unique()
rm(protein_coding_transcripts)
# Location of mutated sequences
mutated_sequences_path <- config$data_folder # To this path add "/{cancer_type}/mutated_sequences/aa/

#######################################
# PUT TOGETHER MUTATED SEQUENCES OF   #
# CURRENT CHROMOSOME FROM ALL TUMORS  #
#######################################
first <- TRUE
for (tumor in tumor_list) {
    # path for the mutated sequences of this specific tumor
    path <- paste0(mutated_sequences_path,"/", tumor, "/mutated_sequences/aa/chr", chr, ".RData")
    # Start concatenation of all chromosomes
    load(path)
    if (first) {
        mutated_sequences <- result
        first <- FALSE
    }else{
        mutated_sequences <- cbind(mutated_sequences, result)
    }
}
# Only keep transcripts with a valid wt sequence
mutated_sequences <- mutated_sequences[which(rownames(mutated_sequences) %in% rownames(sequences_aa)),]
#################################
#       READ IN HAPLOTYPE_ID    # 
#################################

haplotype_IDs <- read_csv(snakemake@input[["IDs"]], col_names = TRUE)
haplotype_IDs <- haplotype_IDs %>% column_to_rownames(var = 'rowname')
# Only keep IDs of haplotypes of interest
haplotype_IDs <- haplotype_IDs %>% filter(ensembl_transcript_id %in% rownames(mutated_sequences))

#############################
#       BUILD VCF FILE      # 
#############################

# Let's put together the first 9 columns of the vcf file which are:
#
# Chromosome name - CHR
# Position - POS
# Variant ID - ENSTID 
# Reference - ENSTID
# Alternative - HAPLOTYPE_ID
# Quality - set to max quality
# Filter - PASS
# info - leave empty
# Format - GT

# Order transcript_positions for chromosome and position, then add "chr" to 
# the chromosome name
transcript_positions <- transcript_positions %>%    arrange(chromosome_name, transcript_start) %>% 
                                                    mutate(chromosome_name = paste0("chr", chromosome_name))
# Create ALT column
alt_ID_vectors <- haplotype_IDs %>% split(.$ensembl_transcript_id) %>% map(pull, haplotype_ID)
alt_ID_col <- lapply(alt_ID_vectors, function(x){
    x <- sort(x)
    if(!str_detect(x[1], "\\.")){
        x <- x[-1]
    }
    paste(paste0("<",x), collapse = ",")
})
alternatives <- data.frame("ALT" = character())
for (vector in alt_ID_col) {
    alternatives <- rbind(alternatives,data.frame("ALT"=vector))
}
alternatives <- alternatives %>% mutate(REF = str_split_i(str_split_i(str_split_i(ALT,",",1),"\\.",1), "<", 2))

vcf_body <- data.frame("CHROM" = transcript_positions$chromosome_name,
                       "POS" = transcript_positions$transcript_start,
                       "ID" = transcript_positions$ensembl_transcript_id,
                       "REF" = transcript_positions$ensembl_transcript_id) %>% dplyr::rename(`#CHROM` = CHROM)
vcf_body <- vcf_body %>% left_join(alternatives, by = "REF")
vcf_body <- cbind(vcf_body, data.frame("QUAL" = rep(60, length(vcf_body$POS)),
                                       "FILTER" = rep("PASS", length(vcf_body$POS)),
                                       "INFO" = rep(".", length(vcf_body$POS)),
                                       "FORMAT" = rep("GT", length(vcf_body$POS))))
# Add ">" before REF alleles annotation as required in Plink2.0 specifications
vcf_body <- vcf_body %>% mutate(REF = paste0("<", REF))
# Remove from the vcf the lines referring to transcripts with no alternative haplotype
vcf_body <- vcf_body %>% filter(!is.na(ALT))
# Remove variants with more than ALT_limit variants 
vcf_body <- vcf_body %>% filter(sapply(ALT, function(x) length(str_split_1(x, ",")) <= ALT_limit))
# For each transcript, check genotype of each sample
res <- mclapply(vcf_body$ID, FUN = function(transcript){
        # Initialize genotype vector
        gt <- c(transcript)
        # Subset haplotype_ID
        transcript_IDs <- haplotype_IDs %>% filter(ensembl_transcript_id == transcript)
        for (sample in colnames(mutated_sequences)) {
            allele_1 <- transcript_IDs[paste0(str_split_i(mutated_sequences[transcript, sample], "-", 1),transcript),1]
            allele_2 <- transcript_IDs[paste0(str_split_i(mutated_sequences[transcript, sample], "-", 2),transcript),1]
            # Account for missing sequences (could be removed by indels)
            if(is.na(allele_1)){
                allele_1 <- "."
            }else{
                if (str_detect(allele_1, pattern = "\\.")) {
                    allele_1 <- str_split_i(allele_1, "\\.",2)
                }else{
                    allele_1 <- 0
                }
            }
            if(is.na(allele_2)){
                allele_2 <- "."
            }else{
                if (str_detect(allele_2, pattern = "\\.")) {
                    allele_2 <- str_split_i(allele_2, "\\.", 2)
                }else{
                    allele_2 <- 0
                }
            }
            gt <- c(gt, paste0(allele_1, "/", allele_2))
        }
        return(gt)
    },
    mc.cores = cores,
    mc.preschedule = TRUE,
    mc.cleanup = TRUE)

# Put togheter results from mclapply
genotypes <- data_frame()
for (item in res) {
    genotypes <- rbind(genotypes, item)
}
colnames(genotypes) <- c("ID", colnames(mutated_sequences))
# Put together fixed part of vcf and genotypes
vcf <- vcf_body %>% left_join(genotypes, by = "ID")
# Save vcf file and haplotype IDs
write_tsv(vcf, file = snakemake@output[["vcf"]])
