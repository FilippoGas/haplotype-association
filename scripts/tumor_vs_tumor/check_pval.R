# Filippo Gastaldello - 16/04/25
#
# The results of the common haplotypes associations analysis show too many 
# associations and with too high p-values.

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(parallel)

# SNAKEMAKE INPUTS ----
results_path <- snakemake@input[["plink_results"]]
haplotype_ID_path <- snakemake@input[["haplotype_ID"]]
# SNAKEMAKE PARAMS ----
vcf_location <- snakemake@params[["vcf_location"]]

tumor_type <- str_split_i(basename(results_path), "\\.", 2)
model_type <- str_split_i(basename(results_path), "\\.", 1)

plot_outdir <- paste0(str_split_i(snakemake@output[[1]], "logistic_results", 1),"plots/", tumor_type, "/", model_type, "/")

model_names <- list("hide-covar"="additive",
                    "dominant"="dominant",
                    "recessive"="recessive")

cores_freq = 30

# READ PLINK RESULTS----
plink_results <- read_tsv(file = results_path);
# Only keep results with pvalue<0.0001 and no error code
plink_results <- plink_results %>% filter(P<0.05, ERRCODE==".")

# READ HAPLOTYPE ID ----
# used to retrieve the number of different haplotypes for each transcript, only keep
# transcripts present in the plink results
haplotype_ID <- read_csv(file = haplotype_ID_path) %>% filter(ensembl_transcript_id %in% plink_results$ID)

# READ VCF FILES----
# Used to retrieve information about the genotypes.
# The vcf files are divided per chromosome, let's put them together and only 
# variants present in the plink results
vcf <- data_frame()
for (chromosome in seq(1,22)) {
    vcf <- rbind(vcf,
                 read_tsv(file = paste0(vcf_location,
                                        "/chr",
                                        chromosome,
                                        "_haplotype.vcf"
                                        )
                          ) %>% filter(ID %in% plink_results$ID)
                 )
}

# CORRELATION ANALYSES----

## pvalue/number of haplotypes per transcript----

# The number of haplotype per transcript can be retrieved from the haplotype_ID df
n_haplotypes <- haplotype_ID %>% summarise(haplotypes = n(), .by = ensembl_transcript_id)

p <- plink_results %>%
        left_join(n_haplotypes, by = join_by(ID==ensembl_transcript_id)) %>% 
        select(ID, P, haplotypes) %>%
        ggplot(aes(x = haplotypes, y = P)) +
            geom_point(size = 0.5) +
            stat_cor(method = "spearman", label.y = 0.035, label.x = 200) +
            ggtitle("Number of haplotype vs P-value")
ggsave(p,
       filename = paste0(plot_outdir, "corr_pval_nhap.png"),
       device = "png")


## pvalue/haplotype frequency on population ----

# The frequency of an haplotype in my population can be extracted from the vcf file
n_samples <- ncol(vcf)-9
# prepare df to store frequencies
hap_freq <- plink_results %>%
                select(ID, A1, P) %>%
                mutate(A1 = str_split_i(A1, "<",2)) %>%
                dplyr::rename(ensembl_transcript_id = ID, haplotype_id = A1)
freq = mclapply(hap_freq$haplotype_id, 
                function(id){ # given a haplotype ID return its frequency
                        # haplotype_id is in the form "<ENSTID.HAP_NUMBER"
                        enst_id <- str_sub(id, 1, 15)
                        hap_number <- ifelse(nchar(id)>15, str_split_i(id, "\\.", 2), "0")
                        # retrieve the vcf row relative to the haplotype's transcript
                        vcf_row <- vcf %>%  filter(ID==enst_id) %>%
                                            select(-c("#CHROM", "POS","ID","REF","ALT","QUAL","FILTER", "INFO","FORMAT")) %>%
                                            as.character() %>%
                                            paste0(collapse = " ")
                        # count how many times the haplotype number occurs in the genotypes
                        count <- str_count(vcf_row, hap_number)
                        return((100*count)/(2*n_samples))
                },
                mc.cores=cores_freq,
                mc.preschedule = TRUE)
# extract list's elements and use them as column for correlation
hap_freq[,"freq"] <- sapply(freq, function(x){ x[[1]] })
p <- hap_freq %>% ggplot(aes(x = freq, y = P)) +
                    geom_point(size = 0.5) +
                    stat_cor(method = "spearman", label.y = 0.045, label.x = 200) +
                    ggtitle("Haplotype frequency vs P-value")
ggsave(p,
       filename = paste0(plot_outdir, "corr_pval_hap_freq.png"),
       device = "png")


# WRITE FAKE SNAKEMAKE OUTPUT
writeLines("done", snakemake@output[[1]])