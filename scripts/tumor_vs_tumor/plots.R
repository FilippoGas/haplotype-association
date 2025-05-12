# Filippo Gastaldello - 15/04/25
#
# Plot results of plink analysis

library(ggplot2)
library(tidyverse)
library(qqman)

# Get input location and read file
plink_result <- snakemake@input[[1]]
plink_result <- read_tsv(plink_result)
logistic_result <- snakemake@input[[2]]
logistic_result <- read_tsv(logistic_result) %>% mutate(ID=str_sub(haplotype_ID,2,16)) %>% dplyr::rename(P="P(>|z|)") %>% transform(P=as.numeric(P))


# Get output path
outdir <- str_split_i(snakemake@output[[1]], ".done", 1)

# Read in transcripts annotations, cancer type and model
anno <- load(snakemake@params[["anno"]])
tumor <- str_split_i(str_split_i(snakemake@output[[1]],"plots/",2),"/done",1)
model <- str_split_i(tumor, "/",2)
tumor <- str_split_i(tumor, "/",1)


n_transcript <- length(unique(plink_result$ID))
plink_result <- plink_result %>% mutate(Bonferroni = P*n_transcript) %>% 
                                dplyr::rename(CHR=`#CHROM`)
# Add affected gene symbol to each entry
plink_result <- left_join(plink_result,
                          protein_coding_transcripts %>% select(ensembl_transcript_id, ensembl_gene_id, hgnc_symbol) %>% unique(),
                          by = join_by(ID==ensembl_transcript_id))
logistic_result <- left_join(logistic_result,
                          protein_coding_transcripts %>% select(ensembl_transcript_id, ensembl_gene_id, hgnc_symbol) %>% unique(),
                          by = join_by(ID==ensembl_transcript_id))
# get chromosome and position for logistic results
logistic_result <- logistic_result %>% left_join(plink_result %>% select(CHR, POS, ID) %>% unique(), by = "ID")

# Annotated qqman manhattan
pdf(file = paste0(outdir, "manhattan_annotated_plink.pdf"), width = 15, height = 10)
manhattan(plink_result %>% 
              filter(P<0.05,
                     ERRCODE=="."),
          chr = "CHR",
          bp = "POS",
          p = "P",
          snp = "hgnc_symbol",
          annotatePval = 0.001,
          genomewideline = -log10(1e-06),
          suggestiveline = -log10(1e-04))
dev.off()

pdf(file = paste0(outdir, "manhattan_annotated_logistic.pdf"), width = 15, height = 10)
manhattan(logistic_result %>% 
              filter(P<0.05,
                     !is.na(P)),
          chr = "CHR",
          bp = "POS",
          p = "P",
          snp = "hgnc_symbol",
          logp = TRUE,
          annotatePval = 0.0001,
          genomewideline = -log10(1e-06),
          suggestiveline = -log10(1e-04))
dev.off()

# Q-Q plot
pdf(file = paste0(outdir, "QQ_plink.pdf"), width = 15, height = 15)
qq(plink_result$P, main = paste0(tumor, "-", model))
dev.off()

pdf(file = paste0(outdir, "QQ_logistic.pdf"), width = 15, height = 15)
qq(logistic_result$P, main = paste0(tumor, "-", model))
dev.off()

# Write fake output
write_file("done", file = snakemake@output[[1]])
