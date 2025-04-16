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

