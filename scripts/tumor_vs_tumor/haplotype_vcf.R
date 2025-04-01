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

################################
# READ INPUTS FROM CONFIG.YAML #
################################
config <- yaml::read_yaml(paste0(snakefile_abs_path,'/config.yaml'))
# Cores for parallel computing (WARNING: vcf creation fase uses 2x cores)
cores = config$cores_make_vcf_tt
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
                            filter(ensembl_transcript_id %in% rownames(sequences_aa)) %>%
                            select(ensembl_transcript_id, transcript_start, chromosome_name) %>% 
                            unique()
rm(protein_coding_transcripts)
# Location of mutated sequences
mutated_sequences_path <- config$data_folder # To this path add "/{cancer_type}/mutated_sequences/aa/

#################################
# CREATE HAPLOTYPE ID DATAFRAME #
#################################

# Each tumor type has the same set of transcripts analysed, but the haplotypes
# may be different from tumor to tumor. For this reason the different files
# containing the different haplotypes for each transcript in each tumor need to 
# be merged together.

# Create df structure for the new IDs
haplotype_IDs <- data.frame("haplotype_ID" = character(),
                            "sequence" = character())

res <- mclapply(rownames(sequences_aa),
                FUN = function(transcript){
                        # Given a transcript look for its file in each tumor's folder and merge 
                        # them together
                        transcript_haplotypes <- data.frame()
                        for (tumor in tumor_list) {
                            transcript_haplotypes <- rbind(transcript_haplotypes,
                                                           read_csv(paste0(config$data_folder, "/",tumor, "/ESM_inputs/",transcript, ".csv"),
                                                                    col_names = FALSE))
                        }
                        # Some tumors might have common haplotypes
                        transcript_haplotypes <- transcript_haplotypes %>% unique() %>% dplyr::rename(haplotypes = X1)
                        # Now that all the unique haplotypes for the current transcripts are collected,
                        # They must be associated to an ID (ENSTID+counter)
                        hap_counter <- 1
                        for (haplotype in transcript_haplotypes$haplotypes) {
                            # The wt sequence will not have a suffix number to be recognized as the 
                            # reference sequence
                            if (!is.na(sequences_aa[transcript,'sequence']) & (haplotype == str_split_i(sequences_aa[transcript, 'sequence'], "\\*", 1))) {
                                haplotype_IDs <- rbind(haplotype_IDs,
                                                       c("haplotype_ID" = transcript,
                                                         "sequence" = haplotype))
                            }else{
                                haplotype_IDs <- rbind(haplotype_IDs,
                                                       c("haplotype_ID" = paste0(transcript, ".", hap_counter),
                                                         "sequence" = haplotype))
                                hap_counter <- hap_counter + 1
                            }
                        }
                        # Reset column names as they get changed
                        colnames(haplotype_IDs) <- c("haplotype_ID", "sequence")
                        return(haplotype_IDs)
                    },
                mc.preschedule = TRUE,
                mc.cores = cores,
                mc.cleanup = TRUE)

# Concatenate results together
haplotype_IDs <- bind_rows(res)

#######################################
# PUT TOGETHER MUTATED SEQUENCES FROM # 
# ALL SAMPLES FROM ALL TUMORS         #
#######################################

res <- mclapply(tumor_list,
                FUN = function(tumor){
                    # path for the mutated sequences of this specific tumor
                    path <- paste0(mutated_sequences_path,"/", tumor, "/mutated_sequences/aa/")
                    # Start concatenation of all chromosomes
                    first <- TRUE
                    for (chr in list.files(path)) {
                        load(paste0(path, chr))
                        if (first) {
                            mutated_sequences <- result
                            first <- FALSE
                        }else{
                            mutated_sequences <- rbind(mutated_sequences, result)
                        }
                    }
                    return(mutated_sequences)
                },
                mc.cores = min(cores, length(tumor_list)),
                mc.preschedule = TRUE,
                mc.cleanup = TRUE)
# Concatenate all tumors together
mutated_sequences <- bind_cols(res)
# Only keep transcripts with a valid wt sequence
mutated_sequences <- mutated_sequences[which(rownames(mutated_sequences) %in% rownames(sequences_aa)),]


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
haplotype_IDs <- haplotype_IDs %>% mutate(ensembl_transcript_id = str_split_i(haplotype_ID, "\\.",1))

# Create ALT column
alt_ID_vectors <- haplotype_IDs %>% split(.$ensembl_transcript_id) %>% map(pull, haplotype_ID)
alt_ID_col <- lapply(alt_ID_vectors, function(x){
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
# Set rownames for haplotype_ID for faster retrieval during the generation of the genotypes
haplotype_IDs <- haplotype_IDs %>% mutate(rowname = paste0(sequence, ensembl_transcript_id))
haplotype_IDs <- haplotype_IDs %>% column_to_rownames(var = "rowname")
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
write_csv(haplotype_IDs, file = snakemake@output[["IDs"]])
