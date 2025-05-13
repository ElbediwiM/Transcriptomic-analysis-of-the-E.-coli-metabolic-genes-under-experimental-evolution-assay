## Download the transcriptomic data for the evolved isolates either from supplementary data or from Gene Expression Omnibus (GEO).

### Filter the transcriptomic data for the metabolic genes using python 

import pandas as pd

# Step 1: Load the CSV file
csv_file_path = 'e:/metabolic_tr/transcriptomic/E.coli/file.csv'  
df = pd.read_csv(csv_file_path)


# Step 2: Identify duplicates in the 'genes' column 
genes_df = pd.read_csv('e:/metabolic_tr/transcriptomic/ E.coli /metabolic_genes.csv')
genes_list = genes_df['gene'].tolist()  

# Step 3: Filter the DataFrame

filtered_df = df[df['gene'].isin(genes_list)]  

# Step 4: Save the filtered DataFrame to a new CSV
output_file_path = 'e:/metabolic_tr/transcriptomic/ E.coli /Filtered data.csv'  
filtered_df.to_csv(output_file_path, index=False)

print(f"Filtered data saved to {output_file_path}")


##### Merge the filtered data of the different studies using R 
# Load necessary libraries
library(tidyverse)

# Step 0: Read the data
data <- read.csv("filtered_genes.csv", header = TRUE)

# Reshape the data into long format
long_data <- data %>%
  pivot_longer(cols = -Genes, names_to = "Condition", values_to = "log2FC")

# Calculate mean log2FC for each gene under each condition
mean_log2FC <- long_data %>%
  group_by(Genes, Condition) %>%
  summarise(mean_log2FC = mean(log2FC, na.rm = TRUE), .groups = 'drop')

# Reshape back to wide format
wide_data <- mean_log2FC %>%
  pivot_wider(names_from = Condition, values_from = mean_log2FC) %>%
  column_to_rownames(var = "Genes")

# Convert to matrix
expression_matrix <- as.matrix(wide_data)

# Replace missing or infinite values with row means
expression_matrix[is.na(expression_matrix)] <- rowMeans(expression_matrix, na.rm = TRUE)[row(expression_matrix)[is.na(expression_matrix)]]
expression_matrix[is.infinite(expression_matrix)] <- rowMeans(expression_matrix, na.rm = TRUE)[row(expression_matrix)[is.infinite(expression_matrix)]]

# Step 1: Define treatment-control pairs for each antibiotic
treatment_control_pairs <- list(
  KAN = c("KAN4", "KAN1"),
  AMP = c("AMP9", "AMP1"),
  CHL = c("CHL9", "CHL1"),
  CPR = c("CPR9", "CPR1"),
  DOX = c("DOX5", "DOX1"),
  ERY = c("ERY6", "ERY1"),
  FOX = c("FOX8", "FOX1"),
  TRM = c("TRM7", "TRM1"),
  NAL = c("NAL8", "NAL1"),
  NIT = c("NIT7", "NIT1"),
  TET = c("TET7", "TET1"),
  TOB = c("TOB9", "TOB1")
)

# Ensure all conditions exist in the data
all_conditions <- unlist(treatment_control_pairs)
if (!all(all_conditions %in% colnames(expression_matrix))) {
  stop("One or more conditions not found in the data.")
}

# Step 2: Compute log2FC, p-values, and q-values for each gene per antibiotic
results_list <- list()

for (antibiotic in names(treatment_control_pairs)) {
  treatment <- treatment_control_pairs[[antibiotic]][1]  # Treatment condition
  control <- treatment_control_pairs[[antibiotic]][2]    # Control condition
  
  # Extract treatment and control values
  if (!(treatment %in% colnames(expression_matrix)) || !(control %in% colnames(expression_matrix))) {
    warning(paste("Skipping", antibiotic, ": Treatment or control condition not found."))
    next
  }
  
  treatment_values <- expression_matrix[, treatment]
  control_values <- expression_matrix[, control]
  
  # Calculate log2 fold change for the current antibiotic
  log2fc <- treatment_values - control_values
  
  # Assume a fixed variance (you can adjust this value based on your data's variability)
  fixed_variance <- 0.1  # Adjust this value as needed
  t_stat <- log2fc / sqrt(fixed_variance)
  
  # Compute p-values using the t-distribution
  pvalue <- 2 * pt(-abs(t_stat), df = ncol(expression_matrix) - 1)
  
  # Adjust p-values for multiple testing using the Benjamini-Hochberg method
  qvalue <- p.adjust(pvalue, method = "BH")
  
  # Store results for the current antibiotic
  results_list[[paste0(antibiotic, "_log2FC")]] <- log2fc
  results_list[[paste0(antibiotic, "_pvalue")]] <- pvalue
  results_list[[paste0(antibiotic, "_qvalue")]] <- qvalue
}

# Step 3: Combine all results into a single data frame
final_results <- do.call(cbind, results_list)


# Add gene names as row names
rownames(final_results) <- rownames(expression_matrix)

# Convert to data frame
final_results <- as.data.frame(final_results)

# Step 4: Reorder columns for clarity
final_results <- final_results %>%
  select(
    starts_with("KAN"), 
    starts_with("AMP"), 
    starts_with("CHL"), 
    starts_with("CPR"), 
    starts_with("DOX"), 
    starts_with("ERY"), 
    starts_with("FOX"), 
    starts_with("TRM"), 
    starts_with("NAL"), 
    starts_with("NIT"), 
    starts_with("TET"), 
    starts_with("TOB")
  )

# Step 5: View top results
head(final_results)

# Step 6: Save results to a CSV file
write.csv(final_results, "gene_expression_per_antibiotic_columns.csv", row.names = TRUE)
