# Filippo Gastaldello - 29/04/2025
#
# Run logistic regression on all significant plink associations to validate them.

library(tidyverse)
library(caret)
library(parallel)

# SNAKEMAKE INPUT ----
results_path <- snakemake@input[["plink_results"]]
haplotype_ID_path <- snakemake@input[["haplotype_ID"]]
covariates_path <- snakemake@input[["covariates"]]
# SNAKEMAKE PARAMS ----
vcf_location <- snakemake@params[["vcf_location"]]
clinical_data_path <- snakemake@params[["clinical_data"]]

tumor_type <- str_split_i(basename(results_path), "\\.", 2)
model_type <- str_split_i(basename(results_path), "\\.", 1)

logistic_outdir <- str_split_i(snakemake@output[[1]], ".done", 1)

model_names <- list("hide-covar"="additive",
                    "dominant"="dominant",
                    "recessive"="recessive")

cores_logistic = 30

# READ PLINK RESULTS----
plink_results <- read_tsv(file = results_path);
# Only keep results with pvalue<0.0001 and no error code
plink_results <- plink_results %>% filter(P<0.05, ERRCODE==".")


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

# LOGISTIC REGRESSION ----

# Run logistic regression to predict cancer type vs all others using significant associations from plink

# Open clinical data and covariates table before mclapply to avoid opening it in all forked instances
clinical_data <- read_tsv(clinical_data_path) %>% select(`Patient ID`, `TCGA PanCanAtlas Cancer Type Acronym`) %>% unique()
covariates <- read_tsv(covariates_path)

logreg_results <- mclapply(plink_results$A1,       # For each significative association
                           function(hap_id){
                               
                               # Extract ENSTID and haplotype number 
                               enst_id <- str_sub(hap_id, 2, 16)
                               hap_number <- ifelse(nchar(hap_id)>15, str_split_i(hap_id, "\\.", 2), "0")
                               
                               # For each haplotype in the plink results, make a table
                               # with sample-condition-genotype-covariates to use in
                               # logistic regression
                               data <- vcf %>% dplyr::filter(ID == enst_id) %>% 
                                   dplyr::select(-c(1:9)) %>%
                                   t() %>% 
                                   as.data.frame() %>% 
                                   rownames_to_column(var = "ID") %>% 
                                   dplyr::rename("GT" = "V1")
                               
                               # Based on the haplotype I am considering change the genotype 
                               # wrt to it in the 3 association model (additive, recessive and dominant)
                               data <- data %>%    dplyr::mutate(additive = str_count(GT, paste0("(?<!\\d)", hap_number,"(?!\\d)")),
                                                                 recessive = ifelse(str_count(GT, paste0("(?<!\\d)", hap_number,"(?!\\d)"))==2, 1, 0),
                                                                 dominant = ifelse(str_count(GT, paste0("(?<!\\d)", hap_number,"(?!\\d)"))==0, 0, 1))
                               
                               # Add condition information from clinical data coded as 
                               # Tumor : 0
                               # All others : 1
                               # and remove samples with NAs
                               data <- data %>%    dplyr::mutate(`Patient ID` = str_sub(ID, 1, 12)) %>%
                                   dplyr::left_join(clinical_data) %>%
                                   dplyr::select(-c(`Patient ID`)) %>% 
                                   dplyr::rename(condition = `TCGA PanCanAtlas Cancer Type Acronym`) %>% 
                                   dplyr::mutate(condition = ifelse(condition==tumor_type, tumor_type, "others")) %>% 
                                   drop_na()
                               # Add covariates
                               data <- data %>% left_join(covariates, by = join_by("ID"=="#IID")) %>% select(-ID) %>% drop_na()
                               
                               # Set condition as factor
                               data$condition <- as.factor(data$condition)
                               current_levels <- levels(data$condition)
                               new_levels <- make.names(current_levels)
                               levels(data$condition) <- new_levels
                               
                               # Split data in train and test (70/30)
                               train_index <- createDataPartition(y = data$condition,
                                                                  p = 0.7,
                                                                  list = FALSE,
                               )
                               train_data <- data[train_index,]
                               test_data <- data[-train_index,]
                               
                               # Specify repeated cross-validation as resampling method
                               ctrl <- trainControl(method = "repeatedcv",
                                                    repeats = 3,
                                                    classProbs = TRUE,
                                                    summaryFunction = twoClassSummary)
                               # Train model
                               trained_model <- train(
                                   as.formula(paste0("condition ~ ", as.name(model_names[[model_type]]), " + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")),
                                   data = data,
                                   method = "glm",
                                   family = "binomial",
                                   preProcess = c("center", "scale"),
                                   trControl = ctrl,
                                   metric = "ROC"
                               )
                               
                               # # Save model summary and results
                               # train_results <- as.data.frame(trained_model$results)
                               # write_tsv(train_results, file = paste0(logistic_outdir, "trained_model_results.tsv"))
                               # train_final_model <- capture.output(trained_model$finalModel)
                               # write_lines(train_final_model, file = paste0(logistic_outdir, "trained_final_model.txt"))
                               # train_model_summary <- capture.output(summary(trained_model))
                               # write_lines(train_model_summary, file = paste0(logistic_outdir, "trained_model_summary.txt"))
                               # # Save model itself
                               # saveRDS(trained_model, file = paste0(logistic_outdir, "trained_model.rds"))
                               # 
                               # # Predict test set
                               # predicted_classes <- predict(trained_model, test_data)
                               # pred_confusion_matrix <- capture.output(confusionMatrix(predicted_classes, test_data$condition))
                               # write_lines(pred_confusion_matrix, file = paste0(logistic_outdir, "pred_conf_matrix.txt"))
                               
                               summary <-  cbind(hap_id,
                                                 as.data.frame(t(str_split(capture.output(summary(trained_model))[12], "\\s+")[[1]])))
                               return(summary[,-length(colnames(summary))])
                               
                           },
                           mc.preschedule = TRUE,
                           mc.cores = cores_logistic)

# Put together results from different models in one dataframe
res <- bind_rows(logreg_results)
colnames(res) <- c("haplotype_ID", "model", "coeff", "std_err", "zvalue", "P(>|z|)", "significance")

# Write snakemake output
write_tsv(res, file = snakemake@output[[1]])