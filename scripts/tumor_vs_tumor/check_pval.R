# Filippo Gastaldello - 16/04/25
#
# The results of the common haplotypes associations analysis show too many 
# associations and with too high p-values.

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(parallel)
library(hrbrthemes)
library(viridis)
library(patchwork)

# SNAKEMAKE INPUTS ----
plink_results_path <- snakemake@input[["plink_results"]]
logistic_results_path <- snakemake@input[["logistic_results"]]
haplotype_ID_path <- snakemake@input[["haplotype_ID"]]
# SNAKEMAKE PARAMS ----
vcf_location <- snakemake@params[["vcf_location"]]

tumor_type <- str_split_i(basename(plink_results_path), "\\.", 2)
model_type <- str_split_i(basename(plink_results_path), "\\.", 1)

plot_outdir <- paste0(str_split_i(snakemake@output[[1]], "plots", 1),"plots/", tumor_type, "/", model_type, "/")

model_names <- list("hide-covar"="additive",
                    "dominant"="dominant",
                    "recessive"="recessive")

cores_freq = snakemake@params[["cores"]]

# READ PLINK AND LOGISTIC REGRESSION RESULTS----
plink_results <- read_tsv(file = plink_results_path)
logreg_results <- read_tsv(file = logistic_results_path) %>%
                    dplyr::rename(P="P(>|z|)") %>% 
                    mutate(ID = str_sub(haplotype_ID, 2 , 16),
                           P = ifelse(P=="<", 2e-16, P),
                           haplotype_id = str_split_i(haplotype_ID, "<", 2)) %>% 
                    transform(P=as.numeric(P)) %>% 
                    select(-haplotype_ID)
# Only keep results with pvalue<0.05 and no error code
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

p1 <- plink_results %>%
        left_join(n_haplotypes, by = join_by(ID==ensembl_transcript_id)) %>% 
        select(ID, P, haplotypes) %>%
        ggplot(aes(x = haplotypes, y = P)) +
            geom_point(size = 0.3) +
            stat_cor(method = "spearman") +
            ggtitle(paste0("Number of haplotype vs P-value (Plink), n = ",nrow(plink_results)))
p2 <- logreg_results %>%
        filter(P<0.05) %>% 
        left_join(n_haplotypes, by = join_by(ID==ensembl_transcript_id)) %>% 
        select(ID, P, haplotypes) %>%
        ggplot(aes(x = haplotypes, y = P)) +
            geom_point(size = 0.3) +
            stat_cor(method = "spearman") +
            ggtitle(paste0("Number of haplotype vs P-value (Logistic), n = ", nrow(logreg_results %>% filter(P<0.05))))
p <- p1 / p2
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
# Merge with results from logistic results
hap_freq  <- hap_freq %>%
                left_join(logreg_results, by = "haplotype_id")%>% 
                dplyr::rename(P_plink=P.x,
                              P_logreg=P.y)
p1 <- hap_freq %>% ggplot(aes(x = freq, y = P_plink)) +
                    geom_point(size = 0.3) +
                    stat_cor(method = "spearman") +
                    ggtitle(paste0("Haplotype frequency vs P-value (Plink), n = ", nrow(hap_freq)))
p2 <- hap_freq %>% ggplot(aes(x = freq, y = P_logreg)) +
                    geom_point(size = 0.3) +
                    stat_cor(method = "spearman") +
                    ggtitle(paste0("Haplotype frequency vs P-value (Logistic), n = ", nrow(hap_freq %>% filter(!is.na(P_logreg)))))
p <- p1 / p2
ggsave(p,
       filename = paste0(plot_outdir, "corr_pval_hap_freq.png"),
       device = "png")


## Plink p-values vs logistic regression p-values ----

pvalues <- plink_results %>%    inner_join(logreg_results %>%
                                            mutate(haplotype_id = paste0("<", haplotype_id)) %>% 
                                            filter(P<0.05), 
                                        by = join_by(A1==haplotype_id)) %>% 
                                select(P.x, P.y) %>% 
                                dplyr::rename(P_plink="P.x",
                                              P_logreg="P.y")
p <- pvalues %>% ggplot(aes(x = P_plink, y = P_logreg)) +
                    geom_point(size=0.3) +
                    stat_cor(method = "spearman") +
                    ggtitle("Correlation between Plink's and logistic regression's p-values (showing P<0.05)")
                    
# P-VAULE DISTRIBUTION----
# Compute Pvalues distribution over AF, at different thresholds of pvalues
AF_breaks <- list(0.5, 0.25, 0.1, 0.05, 0.01, 0.005, 0.001)
pval_thresh <- list(0.05, 1e-03, 1e-04, 1e-05, 1e-06)
pval_dist_plink <- data_frame("AF" = numeric(), "pval_thresh" = character(), "count" = double())
pval_dist_logreg <- data_frame("AF" = numeric(), "pval_thresh" = character(), "count" = double())
for (AF_value in AF_breaks) {
    for (pval in pval_thresh) {
         count_plink <- hap_freq %>% filter(P_plink<pval, freq > AF_value) %>% nrow()
         count_logreg <- hap_freq %>% filter(P_logreg<pval, freq > AF_value) %>% nrow()
         pval_dist_plink <- rbind(pval_dist_plink, c(AF_value, as.character(pval), as.double(count_plink)))
         pval_dist_logreg <- rbind(pval_dist_logreg, c(AF_value, as.character(pval), as.double(count_logreg)))
    }
}
colnames(pval_dist_plink) <- c("AF", "pval_thresh", "count")
colnames(pval_dist_logreg) <- c("AF", "pval_thresh", "count")

p1 <- pval_dist_plink %>% ggplot(aes(x = AF, y = as.numeric(count), group = pval_thresh, color = pval_thresh)) +
                            geom_line() +
                            geom_point() +
                            scale_color_viridis(discrete = TRUE) +
                            ggtitle("Haplotype frequency vs P-value (Plink)")
p2 <- pval_dist_logreg %>% ggplot(aes(x = AF, y = as.numeric(count), group = pval_thresh, color = pval_thresh)) +
                            geom_line() +
                            geom_point() +
                            scale_color_viridis(discrete = TRUE) +
                            ggtitle("Haplotype frequency vs P-value (Logistic)")
p <- p1 / p2
ggsave(p,
       filename = paste0(plot_outdir, "pval_dist_over_AF.pdf"),
       device = "pdf")

# WRITE FAKE SNAKEMAKE OUTPUT
writeLines("done", snakemake@output[[1]])
