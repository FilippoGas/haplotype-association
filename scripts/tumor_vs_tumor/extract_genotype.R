# Filippo Gastaldello - 30/04/25
#
# Given a table of mutated sequences and their haplotype ID, extract the genotype
# leading to that sequence

library(tidyverse)
library(vcfR)
library(parallel)

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
results <- c()
for (chr in seq(1,22)) {
    
    # Subset haplotype_ID to the transcripts belonging to this chromosome and remove reference haplotype
    haplotype_ID_chr <- haplotype_ID %>% filter(chromosome_name == chr,
                                                str_detect(haplotype_ID, "\\."))
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
            temp <- result
            first <- FALSE
        }else{
            temp <- cbind(temp, result)
        }
    }
    mutated_sequences <- temp
    rm(temp)
    rm(result)
    
    # Only keep valid mutated sequences
    mutated_sequences_chr <- mutated_sequences[rownames(mutated_sequences) %in% haplotype_ID_chr$ensembl_transcript_id,]
    
    # For each mutated sequence find the first sample with the same sequence and
    # look for his genotype in that specific region
    res <- mclapply(haplotype_ID_chr$haplotype_ID[1:20],
                    function(id) {
                        # Get the transcript ID relative to this haplotype
                        transcript <- haplotype_ID_chr %>% filter(haplotype_ID == id) %>% select(ensembl_transcript_id) %>% as.character()
                        # Get the haplotype sequence
                        hap <- haplotype_ID_chr %>% filter(haplotype_ID == id) %>% select(sequence) %>% as.character()
                        for (sample in colnames(mutated_sequences_chr)) {
                            # Sequences are saved as "seq1-seq2"
                            sequences <- str_split_1(mutated_sequences_chr[transcript, sample], pattern = "-")
                            # Exit loop when a sample with this haplotype is found, and keep track of the allele
                            # containing the desired sequence
                            if (sequences[1]==hap) {
                                allele <- 1
                                break
                            }
                            if (sequences[2]==hap) {
                                allele <- 2
                                break
                            }
                        }
                        # sample found, read vcf from system, regions can be retrieved from
                        # protein_coding_transcripts' exon coordinates
                        exon_regions <- protein_coding_transcripts %>% 
                                            filter(ensembl_transcript_id==transcript)%>% 
                                            mutate(region = paste0("chr", chromosome_name, ":", exon_chrom_start, "-", exon_chrom_end)) %>% 
                                            select(region)
                        exon_regions <- paste0(exon_regions$region, collapse = ",")
                        # Prepare bash command and collect result
                        query <- paste0("bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%GT\n]' -r ",exon_regions, " -s ", sample, " ",  phased_vcf_path, "chr", chr, ".phased.vcf.gz")
                        genotypes <- read.table(text = system(query, intern = TRUE))
                        colnames(genotypes) <- c("CHROM","POS", "REF", "ALT", "GT")
                        # Only keep genotypes, different from REF, for the allele that generate the sequence of interest
                        genotypes <- genotypes %>% mutate(GT = str_split_i(GT, "\\|", as.numeric(allele))) %>% filter(!GT==0)           #CHECK FOR HAPLOTYPES WITH NO VARIANTS
                        # Reformat variants in the fom CHR:POS.REF>ALT
                        genotypes <- genotypes %>% mutate(variants = paste0(CHROM, ":", POS, ".", REF, ">", ALT))
                        
                        # Return variants in vector form
                        return(data_frame("haplotype_ID" = id,
                                          "variants" = paste0(genotypes$variants, collapse = ",")))
                    },
                    mc.cores = 3,
                    mc.preschedule = TRUE,
                    mc.cleanup = TRUE)
    
    res <- bind_rows(res)
    results <- c(results, res)
    
}
