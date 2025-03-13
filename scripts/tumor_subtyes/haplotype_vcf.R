# Filippo Gastaldello - 06/02/2025
#
# Generate a vcf file to describe presence/absence of specific haplotypes
#(mutated sequences) in a population, starting from a sample/sequence matrix.

library(tidyverse)

# A list of all the mutated sequences for each transcript is present in the 
# ESM_input files.
ESM_input_path = snakemake@input[["ESM_input"]]

# Maximum number of alternative alleles allowed per variant. Variants with more 
# ALT alleles will be removed
ALT_limit <- 254

# Get tumor type
tumor_type <- strsplit(basename(snakemake@output[[1]]), "_", 1)
# Import mutated sequences, wt sequences and transcripts annotations.
# Mutated sequences are divided by chromosome so they need to be loaded independently 
# and concatenated

mutated_sequences_path <- snakemake@input[["mutated_sequences"]]
first <- TRUE
for (file in list.files(path = mutated_sequences_path)) {
    
    load(paste0(mutated_sequences_path, "/", file))
    
    if (first) {
        
        mutated_sequences <- result
        first <- FALSE
    
    }else{
        
        mutated_sequences <- rbind(mutated_sequences, result)    
        
    }
    
}

load("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/resources/protein_coding_transcripts/protein_coding_transcripts.RData")
load("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/resources/protein_coding_transcripts/wt_aa.RData")
sequences_aa <- column_to_rownames(sequences_aa, var = "ensembl_transcript_id")

# From the transcript annotation df I only need the transcript start and chromosome
transcript_positions <- protein_coding_transcripts %>% 
    filter(ensembl_transcript_id %in% rownames(mutated_sequences)) %>% 
    dplyr::select(ensembl_transcript_id, transcript_start, chromosome_name) %>% 
    unique()
rm(protein_coding_transcripts)

# Before building any data structure I need new IDs for each haplotype, and a mapping between every 
# haplotype and the respective ID
#
# Create the df structure for the new IDs
haplotypes_ID <- data.frame("haplotype_ID" = character(),
                            "sequence" = character())
colnames(haplotypes_ID) <- c("haplotype_ID", "sequence")

# A unique list of haplotypes for each transcript is present in the ESM_input files
files = list.files(path = ESM_input_path)

# Initialize progress bar to follow progress
pb = txtProgressBar(min = 0, max = length(files), initial = 0, style = 3, title = paste0(tumor_type, " haplotype ID generation"))
i <- 0
for (file in files){
    
    # First, get the ensembl transcript ID from the filename
    enst_id <- str_split_i(file, pattern = ".csv", 1)
    
    # Each file contains all the different haplotypes for one transcript in the population of study
    # cycle them all and add them in the haplotypes_ID df, adding a counter to the original ENSTID
    transcript_haplotypes <- read.csv(paste0(ESM_input_path,"/",file), header = FALSE)
    hap_counter <- 1
    
    for (haplotype in transcript_haplotypes$V1) {
        
        # Check if sequence is wt
        if (!is.na(sequences_aa[enst_id,1]) & (haplotype == str_split_i(sequences_aa[enst_id,1], "\\*",1))) {
            
            haplotypes_ID <- bind_rows(haplotypes_ID,
                                       c("haplotype_ID" = enst_id,
                                         "sequence" = haplotype))
            
        }else{
            
            haplotypes_ID <- bind_rows(haplotypes_ID,
                                       c("haplotype_ID" = paste0(enst_id,".", hap_counter),
                                         "sequence" = haplotype))
            hap_counter <- hap_counter + 1
            
        }
    }
    i <- i + 1
    setTxtProgressBar(pb, i)
}
close(pb)


haplotypes_ID <- haplotypes_ID %>% mutate(ensembl_transcript_id = str_split_i(haplotype_ID, "\\.",1),
                                          rowname = paste0(sequence, ensembl_transcript_id))
haplotypes_ID <- haplotypes_ID %>% column_to_rownames(var = "rowname")

# Add "chr" tp chromosome name in transcript_positions and order for position
transcript_positions <- transcript_positions %>% mutate(chromosome_name = paste0("chr", chromosome_name))
transcript_positions <- transcript_positions %>% arrange(chromosome_name, transcript_start)

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

alt_ID_vectors <- haplotypes_ID %>% split(.$ensembl_transcript_id) %>% map(pull, haplotype_ID)
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

# Create sample genotype section of the vcf
vcf_body <- cbind(vcf_body, as_data_frame(matrix(ncol = length(colnames(mutated_sequences)), nrow = length(rownames(mutated_sequences)))))
colnames(vcf_body)[10:length(colnames(vcf_body))] <- colnames(mutated_sequences) 

# Remove from the vcf the lines referring to transcripts with no alternative haplotype
vcf_body <- vcf_body %>% filter(!is.na(ALT))
# Remove variants with more than ALT_limit variants 
vcf_body <- vcf_body %>% filter(sapply(ALT, function(x) length(str_split_1(x, ",")) <= ALT_limit))

# For each transcript, check genotype of each sample
transcript_count <- 1

# Initiate progress bar to follow progress
pb = txtProgressBar(min = 0, max = length(vcf_body$ID), initial = 0, style = 3, title = paste0(tumor_type, " vcf generation"))

for (transcript in vcf_body$ID) {
    
    
    # Subset haplotype_ID
    transcript_IDs <- haplotypes_ID %>% filter(ensembl_transcript_id == transcript)
    
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
        vcf_body[transcript_count, sample] <- paste0(allele_1, "/", allele_2)
    }
    
    transcript_count <- transcript_count + 1
    
    # Update progress bar
    setTxtProgressBar(pb, transcript_count)
}

close(pb)

# Save vcf as tsv
write_tsv(vcf_body, file = paste0(snakemake@output[["vcf"]]))