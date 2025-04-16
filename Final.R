# Differential Analysis on p53
# Data set is public on kaggle

library(dplyr)
library(ggplot2)
library(DESeq2)
library(tidyverse)

expression_data <- read.csv("C:/Users/tiffa/OneDrive/Desktop/Masters in Bioinformatics/Proteomics/77_cancer_proteomes_CPTAC_itraq.csv/77_cancer_proteomes_CPTAC_itraq.csv", header = TRUE)
head(expression_data) # just to look

# Only want TP53
gene_of_interest <- "NP_000537"
gene_data <- expression_data %>%
  filter(RefSeq_accession_number == gene_of_interest)

head(gene_data) # just to look

# Extract expression values (columns after the first three columns)
expression_data <- gene_data[, 4:ncol(gene_data)]
n_samples <- ncol(expression_data)
# From the data card, the last 3 columns are from healthy individuals
tumor_samples <- expression_data[, 1:(n_samples -3)]
healthy_samples <- expression_data[(n_samples-2): n_samples]

# Convert to numeric vectors
tumor_values <- as.numeric(tumor_samples[1, ])
healthy_values <- as.numeric(healthy_samples[1, ])
# Two-sample t-test (Welch's by default, doesn't assume equal variance)
t_test_result <- t.test(tumor_values, healthy_values)
print(t_test_result)

# plot and reshape
plot_data <- data.frame(
  Expression = c(tumor_values, healthy_values),
  Group = c(rep("Tumor", length(tumor_values)), rep("Healthy", length(healthy_values)))
)

# Violin plot would be best for visualization for single gene
ggplot(plot_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  labs(title = "TP53 Expression: Tumor vs. Healthy",
       y = "Log2 iTRAQ Expression", x = "") +
  theme_minimal() +
  scale_fill_manual(values = c("Tumor" = "#d95f02", "Healthy" = "#1b9e77")) +
  theme(text = element_text(size = 14))

# doing the differential expression for the whole proteome:
expression_data <- read.csv("C:/Users/tiffa/OneDrive/Desktop/Masters in Bioinformatics/Proteomics/77_cancer_proteomes_CPTAC_itraq.csv/77_cancer_proteomes_CPTAC_itraq.csv", header = TRUE)

# Determine the number of columns
n_samples <- ncol(expression_data) - 3  # since 3 last columns are healthy

# Separate tumor and healthy sample columns
tumor_samples <- expression_data[, 4:(ncol(expression_data) - 3)]
healthy_samples <- expression_data[, (ncol(expression_data) - 2):ncol(expression_data)]

# Initialize result dataframe
results <- data.frame(
  gene_symbol = expression_data$gene_symbol,
  log2FC = NA,
  pvalue = NA
)

# Loop through each row (protein)
for (i in 1:nrow(expression_data)) {
  tumor_values <- as.numeric(tumor_samples[i, ])
  healthy_values <- as.numeric(healthy_samples[i, ])
  
  # Skip if too many NAs
  if (sum(!is.na(tumor_values)) >= 3 & sum(!is.na(healthy_values)) >= 3) {
    # t-test
    ttest <- t.test(tumor_values, healthy_values)
    
    # Log2 Fold Change: tumor mean - healthy mean
    log2fc <- mean(tumor_values, na.rm=TRUE) - mean(healthy_values, na.rm=TRUE)
    
    # Store results
    results$log2FC[i] <- log2fc
    results$pvalue[i] <- ttest$p.value
  }
}

# Adjust p-values using FDR)
results$padj <- p.adjust(results$pvalue, method = "BH")

# Add a column to label significant proteins
results$significant <- ifelse(results$padj < 0.05 & abs(results$log2FC) > 1, "Yes", "No")

# Volcano Plot
ggplot(results, aes(x = log2FC, y = -log10(padj))) +
  geom_point(aes(color = significant), alpha = 0.7) +
  scale_color_manual(values = c("No" = "gray", "Yes" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Proteins",
       x = "log2 Fold Change (Tumor vs Healthy)",
       y = "-log10(FDR adjusted p-value)") +
  geom_vline(xintercept = c(-1, 1), linetype="dashed", color="black") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="black")

# Plot looks good, lets look at the top hits
top_hits <- results %>%
  filter(significant == "Yes") %>%
  arrange(padj) 
# Don't see tp53 or cyld in it...
