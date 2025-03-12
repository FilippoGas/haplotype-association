# Filippo Gastaldello - 22/01/2025
#
# Generate a PED/MAP-like file containing haplotypes instead of variants.
# Needs as input a dataframe with samples (patients) on columns and ensembl transcript IDs
# on the rows. The cells will contain the patient-specific sequence for that transcript.
# Create a PED/MAP file to describe the presence/absence of each transcript in each sample.

library(tidyverse)

# ESM_input location
ESM_input_path = "/shares/CIBIO-Storage/BCG/scratch/proteinModel/phasing/ESM_inputs"

output_dir = "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/data/haplotype_PED_MAP/cancer_genes_BRCA/"

# Import mutated sequences and transcripts annotations
load("/shares/CIBIO-Storage/BCG/scratch/proteinModel/phasing/mutated_aa_sequences/mutated_aa_sequences_cancer_genes.RData")
load("/shares/CIBIO-Storage/BCG/scratch/proteinModel/phasing/resources/cds/transcript_annotations.RData")

# From the transcript annotation df I only need the transcript start and chromosome
transcript_positions <- protein_coding_transcripts %>% 
                        filter(ensembl_transcript_id %in% rownames(mutated_sequences)) %>% 
                        select(ensembl_transcript_id, transcript_start, chromosome_name) %>% 
                        unique()
rm(protein_coding_transcripts)

# Before building any data structure I need new IDs for each haplotype, and a mapping between every 
# haplotype and the respective ID
#
# Create the df structure for the new IDs
haplotypes_ID <- data.frame("haplotype_ID" = character(),
                            "sequence" = character(),
                            "start" = character(),
                            "chr" = character())
colnames(haplotypes_ID) <- c("haplotype_ID", "sequence", "start", "chr")


# A unique list of haplotypes for each transcript is present in the ESM_input files
for (file in list.files(path = ESM_input_path)){
    
    # First, get the ensembl transcript ID
    enst_id <- str_split_i(file, pattern = ".csv", 1)
    
    # Each file contains all the different haplotypes for one transcript in the population of study
    # cycle them all and add them in the haplotypes_ID df, adding a counter to the original ENSTID
    transcript_haplotypes <- read.csv(paste0(ESM_input_path,"/",file), header = FALSE)
    hap_counter <- 1
    
    for (haplotype in transcript_haplotypes$V1) {
        
        haplotypes_ID <- bind_rows(haplotypes_ID,
                                   c("haplotype_ID" = paste0(enst_id,".", hap_counter),
                                     "sequence" = haplotype))
        
        hap_counter <- hap_counter + 1
        
    }
}

# Add "chr" tp chromosome name in haplotype_ID and transcript_positions
haplotypes_ID <- haplotypes_ID %>% mutate(chr = paste0("chr", chr))
transcript_positions <- transcript_positions %>% mutate(chromosome_name = paste0("chr", chromosome_name))

# The MAP file is just a subset of transcripts_position's columns in different order
MAP <- data.frame("chromosome" = transcript_positions$chromosome_name,
                  "ID" = transcript_positions$ensembl_transcript_id,
                  "genomic distance" = rep(0, nrow(transcript_positions)),
                  "position" = transcript_positions$transcript_start)

# Extract list of unique ensembl transcript IDs from the MAP files to produce
# a PED file with the same order
transcript_order <- MAP %>% select(ID) %>% mutate(ID = str_split_i(ID, "\\.", 1)) %>% unique() %>% as_vector()


# I can create the backbone of the PED file, composed by the first 6 columns which
# are:
# Family ID (if unknown use the same id as for the sample id in column two)
# Sample ID
# Paternal ID (if unknown use 0)
# Maternal ID (if unknown use 0)
# Sex (if unknown use 0)
# Not used, set to 0

PED <- data.frame("family ID" = colnames(mutated_sequences),
                  "sample ID" = colnames(mutated_sequences),
                  "paternal ID" = rep(0, ncol(mutated_sequences)),
                  "maternal ID" = rep(0, ncol(mutated_sequences)),
                  "sex" = rep(0, ncol(mutated_sequences)),
                  "not used" = rep(0, ncol(mutated_sequences)))

# Cycle transcripts and for each sample check which haplotype is present and fill
# PED file columns accordingly

for (transcript in transcript_order) {

    # In the PED file the two alleles are coded in two adjacent columns
    genotypes_allele_1 <- c()
    genotypes_allele_2 <- c()
    
    for (sample in PED$sample.ID) {
        
        # Each cell in the mutated_sequences dataset are two sequences separated
        # by a "-"
        genotypes_allele_1 <- c(genotypes_allele_1, haplotypes_ID[which((haplotypes_ID$sequence == str_split_i(mutated_sequences[transcript, sample], "-", 1)) & ((str_split_i(haplotypes_ID$haplotype_ID, "\\.", 1) == transcript))),1])
        genotypes_allele_2 <- c(genotypes_allele_2, haplotypes_ID[which((haplotypes_ID$sequence == str_split_i(mutated_sequences[transcript, sample], "-", 2)) & ((str_split_i(haplotypes_ID$haplotype_ID, "\\.", 1) == transcript))),1])
        
    }
    
    # Add transcripts genotypes ar column of the PED file
    
    PED[,paste0(transcript, "_allele_1")] <- genotypes_allele_1
    PED[,paste0(transcript, "_allele_2")] <- genotypes_allele_2
    
}


# Save PED and MAP files