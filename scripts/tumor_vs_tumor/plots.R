# Filippo Gastaldello - 15/04/25
#
# Plot results of plink analysis

library(ggplot2)
library(tidyverse)
library(qqman)

# Get input location and read file
plink_result <- snakemake@input[[1]]
plink_result <- read_tsv(plink_result)

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

# Plot manhattan
don <- plink_result %>% 
        filter(P<0.05,
               ERRCODE==".") %>% 
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(POS)) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(plink_result %>% filter(P<0.05,
                                 ERRCODE=="."), ., by=c("CHR"="CHR")) %>%
    # Add a cumulative position of each SNP
    arrange(CHR, POS) %>%
    mutate( BPcum=POS+tot)

axisdf = don %>%
    group_by(CHR) %>%
    summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

p <- ggplot(don, aes(x=BPcum, y=-log10(P))) +
    
        # Show all points
        geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
        scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
        
        # Add exome-wide and suggested significance level
        geom_hline(yintercept = -log10(1e-06)) +
        geom_hline(yintercept = -log10(1e-04)) +
        annotate("text", x = 0, y=c(-log10(1e-06)+0.5, -log10(1e-04))+0.5, label=c("1e-06", "1e-04")) +
        
        # custom X axis:
        scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
        scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
        
        # Custom the theme:
        theme_bw() +
        theme( 
            legend.position="none",
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()
        ) +
        ggtitle(paste0(tumor, ", ", model, " model - Unadjsuted P-value"))

ggsave(filename = paste0(outdir, "manhattan.pdf"),
       plot = p,
       device = "pdf",
       width = 15,
       height = 10,
       create.dir = TRUE )

# Annotated qqman manhattan
pdf(file = paste0(outdir, "manhattan_annotated.pdf"), width = 15, height = 10)
manhattan(plink_result %>% 
              filter(P<0.05,
                     ERRCODE=="."),
          chr = "CHR",
          bp = "POS",
          p = "P",
          snp = "hgnc_symbol",
          annotatePval = 0.0001,
          genomewideline = -log10(1e-06),
          suggestiveline = -log10(1e-04))
dev.off()

# Q-Q plot
pdf(file = paste0(outdir, "QQ.pdf"), width = 15, height = 15)
qq(don$P, main = paste0(tumor, "-", model))
dev.off()

# Write fake output
write_file("done", file = snakemake@output[[1]])
