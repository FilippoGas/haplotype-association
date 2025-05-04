# Filippo Gastaldello - 30/04/25
#
# Given a table of mutated sequences and their haplotype ID, extract the genotype
# leading to that sequence

library(tidyverse)
library(vcfR)

# SNAKEMAKE INPUT ----
haplotype_ID <- read.csv("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/data/TCGA_haplotype_association/TCGA/tumor_vs_tumor/haplotypes_IDs.csv")
load("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/resources/protein_coding_transcripts/protein_coding_transcripts.RData")

# SNAKEMAKE PARAMS ----
phased_vcf_path <- "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/data/MutSeqGenerator/temp/phased_vcf/"
tumor_types <- c("DLBC","ACC","BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS")

# Let's add chromosome information to haplotype_ID
haplotype_ID <- haplotype_ID %>% left_join(protein_coding_transcripts %>%
                                               select(ensembl_transcript_id, chromosome_name) %>% 
                                               unique(),
                                           by = "ensembl_transcript_id")  

# Let's go through haplotype_ID per chromosome
for (chr in seq(1,22)) {
    
    # Subset haplotype_ID to the transcripts belonging to this chromosome
    haplotype_ID_chr <- haplotype_ID %>% filter(chromosome_name == chr)
    # Put together mutated sequences of this chromosome from all tumors
    mutated_sequences_path <- lapply(tumor_types,
                                     function(x) {
                                         paste0("/shares/CIBIO-Storage/BCG/scratch/fgastaldello/data/MutSeqGenerator/tumors/",
                                                x,
                                                "/mutated_sequences/aa/chr",
                                                chr,
                                                ".RData")
                                     })
    # Open and bind all files
    first <- TRUE
    for (file in mutated_sequences_path) {
        load(file)
        if (first) {
            temp <- mutated_sequences
            first <- FALSE
        }else{
            temp <- cbind(temp, mutated_sequences)
        }
    }
    mutated_sequences <- temp
    rm(temp)
    
    # For each mutated sequence find the first sample with the same sequence and
    # look for his genotype in that specific region
    res <- mclapply(haplotype_ID_chr$sequence,
                    function(hap) {
                        # Get the transcript ID relative to this haplotype
                        transcript <- haplotype_ID %>% filter(sequence == hap) %>% select(ensembl_transcript_id) %>% as.character()
                        for (sample in colnames(mutated_sequences)) {
                            # Sequences are saved as "seq1-seq2"
                            sequences <- str_split_1(mutated_sequences[transcript, sample], pattern = "-")
                            # Exit loop when a sample with this haplotype is found
                            if (sequences[1]==hap) {
                                break
                            }
                            if (sequences[2]==hap) {
                                break
                            }
                        }
                        # sample found, read vcf from system, regions can be retrieved from
                        # protein_coding_transcripts' exon coordinates
                        exon_regions <- protein_coding_transcripts %>% 
                                            filter(ensembl_transcript_id==transcript) %>% 
                                            select(chromosome_name, exon_chrom_start, exon_chrom_end) %>% 
                                            mutate(region = paste0(chromosome_name, ":", exon_chrom_start, "-", exon_chrom_end)) %>% 
                                            select(region)
                        exon_regions <- paste0(exon_regions$region, collapse = ",")
                        # Prepare bash command 
                        query <- paste0("bcftools query ")
                        
                    },
                    mc.cores = 50,
                    mc.preschedule = TRUE,
                    mc.cleanup = TRUE)
    
}
