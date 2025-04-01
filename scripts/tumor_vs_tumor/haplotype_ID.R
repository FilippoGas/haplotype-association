# Filippo gastaldello - 01/04/2025
#
# Generate unique IDs for the haplotype that are going to be genotyped in the vcf
# for the association analysis

library(tidyverse)
library(parallel)

# Get location of snakefile
snakefile_abs_path <- gsub("[\r\n]", "", snakemake@params[["abs_path"]])

################################
# READ INPUTS FROM CONFIG.YAML #
################################
config <- yaml::read_yaml(paste0(snakefile_abs_path,'/config.yaml'))
# Cores for parallel computing
cores = snakemake@params[["cores"]]
# List of tumor types to take into consideration
tumor_list <- config$cancer_types
# wt sequence of proteins of interest, only keep those with a valid sequence
load(paste0(config$wt_transcripts, "/wt_aa.RData"))
sequences_aa <- column_to_rownames(sequences_aa, var = "ensembl_transcript_id")
sequences_aa <- sequences_aa %>% filter(sequence != "Sequence unavailable")

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

haplotype_IDs <- haplotype_IDs %>% mutate(ensembl_transcript_id = str_split_i(haplotype_ID, "\\.",1))
# Set rownames for haplotype_ID for faster retrieval during the generation of the genotypes
haplotype_IDs <- haplotype_IDs %>% mutate(rowname = paste0(sequence, ensembl_transcript_id))

write_csv(haplotype_IDs, file = snakemake@output[["IDs"]])