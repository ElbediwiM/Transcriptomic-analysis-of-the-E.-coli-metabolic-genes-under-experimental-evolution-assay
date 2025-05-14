##### Compute log2FC, P- value and Q value of the filtered data the using R 
# load the libraries

library(statmod)
library(limma)


# Load your count matrix and metadata
count_data <- read.csv("filter_genes1.csv", row.names = 1)
metadata <- read.csv("meta_data.csv", row.names = 1)


# Read the metadata CSV file into a data frame
metadata <- read.csv("meta_data.csv", header = TRUE, stringsAsFactors = TRUE)


metadata$gene <-factor(metadata$gene)
metadata$Condition <- factor(metadata$Condition)


design <- model.matrix(~ gene + Condition, data = metadata)


# Ensure the column names of the count matrix match the row names of the metadata
all(colnames(count_data) == rownames(metadata))  # Should return TRUE


design <- model.matrix(~ gene * Condition, data = metadata)



metadata$Condition <- factor(make.names(metadata$Condition))
levels(metadata$Condition)

design <- model.matrix(~0 + Condition, data = metadata)
colnames(design) <- levels(metadata$Condition)  # Ensure valid names

contrast.matrix <- makeContrasts(
  AMX = Treatment_AMX - Control_AMX,
  CEF = Treatment_CEF - Control_CEF,
  CFT = Treatment_CFT - Control_CFT,
  CIP = Treatment_CIP - Control_CIP,
  COT = Treatment_COT - Control_COT,
  DAP = Treatment_DAP - Control_DAP,
  IMI = Treatment_IMI - Control_IMI,
  KAN = Treatment_KAN - Control_KAN,
  LIN = Treatment_LIN - Control_LIN,
  LVX = Treatment_LVX - Control_LVX,
  MOX = Treatment_MOX - Control_MOX,
  PEN = Treatment_PEN - Control_PEN,
  RIF = Treatment_RIF - Control_RIF,
  TET = Treatment_TET - Control_TET,
  TOB = Treatment_TOB - Control_TOB,
  VNC = Treatment_VNC - Control_VNC,
  levels = design
)


fit <- lmFit(count_data, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results <- topTable(fit2, coef="AMX", adjust="fdr", number=Inf)  # Example for AMX


# Create a list to store results
results_list <- list()

# Extract results for each antibiotic and store them
for (antibiotic in colnames(contrast.matrix)) {
  results_list[[antibiotic]] <- topTable(fit2, coef = antibiotic, adjust = "fdr", number = Inf)
  
  # Save each result as a CSV file
  write.csv(results_list[[antibiotic]], file = paste0("DEG_results_", antibiotic, ".csv"))
}

# Print summary
print("Differential expression analysis completed. Results saved as CSV files.")

# Combine all results into one data frame
all_results_filtered_genes1 <- do.call(rbind, lapply(names(results_list), function(antibiotic) {
  df <- results_list[[antibiotic]]
  df$Antibiotic <- antibiotic  # Add antibiotic as a new column
  return(df)
}))

# Ensure logFC and p-value columns are present
all_results_filtered_genes1$logP <- -log10(all_results$P-Value)


###### Merge log2FC, P- value and Q value of the filtered data of the different studies using R
# Load necessary libraries
library(dplyr)
library(tidyr)

# Load the first CSV file
df1 <- read.csv("all_results_filtered_genes1.csv")

# Load the second CSV file
df2 <- read.csv("all_results_filtered_genes2.csv")


# Melt the first dataframe to long format
df1_melted <- df1 %>%
  pivot_longer(cols = starts_with("AMP") | starts_with("TET") | starts_with("KAN") | 
                 starts_with("CHL") | starts_with("CPR") | starts_with("DOX") | 
                 starts_with("ERY") | starts_with("FOX") | starts_with("TRM") | 
                 starts_with("NAL") | starts_with("NIT") | starts_with("TOB"),
               names_to = "Metric",
               values_to = "Value")

# Split the Metric column into Antibiotic and Metric Type
df1_melted <- df1_melted %>%
  separate(Metric, into = c("Antibiotic", "Metric_Type"), sep = "_")

# Pivot the melted dataframe to wide format
df1_pivot <- df1_melted %>%
  pivot_wider(names_from = Metric_Type, values_from = Value)

# Melt the second dataframe to long format
df2_melted <- df2 %>%
  pivot_longer(cols = -Genes,
               names_to = "Metric",
               values_to = "Value")

# Split the Metric column into Antibiotic and Metric Type
df2_melted <- df2_melted %>%
  separate(Metric, into = c("Antibiotic", "Metric_Type"), sep = "_")

# Pivot the melted dataframe to wide format
df2_pivot <- df2_melted %>%
  pivot_wider(names_from = Metric_Type, values_from = Value)

# Merge the two dataframes on Genes and Antibiotic
combined_df <- full_join(df1_pivot, df2_pivot, by = c("Genes", "Antibiotic"), suffix = c("_file1", "_file2"))

# Function to calculate weighted average for log2FC
weighted_average_log2fc <- function(log2fc1, log2fc2, pvalue1, pvalue2) {
  # Handle NA values
  if (is.na(log2fc1)) log2fc1 <- 0
  if (is.na(log2fc2)) log2fc2 <- 0
  if (is.na(pvalue1)) pvalue1 <- Inf  # Set to Inf to give it zero weight
  if (is.na(pvalue2)) pvalue2 <- Inf  # Set to Inf to give it zero weight
  
  # Use inverse variance as weights (1 / p-value)
  weight1 <- ifelse(pvalue1 > 0, 1 / pvalue1, 0)
  weight2 <- ifelse(pvalue2 > 0, 1 / pvalue2, 0)
  
  # Calculate weighted average
  if (weight1 + weight2 > 0) {
    return((log2fc1 * weight1 + log2fc2 * weight2) / (weight1 + weight2))
  } else {
    return(NA)
  }
}

# Function to calculate geometric mean
geometric_mean <- function(values) {
  values <- values[!is.na(values) & values > 0]  # Ignore NAs and zeros
  if (length(values) == 0) {
    return(NA)
  }
  return(exp(mean(log(values))))
}

# Apply the functions to calculate combined values
combined_df <- combined_df %>%
  rowwise() %>%
  mutate(
    log2FC_combined = weighted_average_log2fc(log2FC_Filtered data1, log2FC_Filtered data2, pvalue_Filtered data1, pvalue_Filtered data2),
    pvalue_combined = geometric_mean(c(pvalue_Filtered data1, pvalue_Filtered data2)),
    qvalue_combined = geometric_mean(c(qvalue_Filtered data1, qvalue_Filtered data2))
  ) %>%
  ungroup()

# Select relevant columns
combined_df <- combined_df %>%
  select(Genes, Antibiotic, log2FC_combined, pvalue_combined, qvalue_combined)

# Pivot the combined data to wide format with antibiotics in columns
combined_wide <- combined_df %>%
  pivot_wider(names_from = Antibiotic, 
              values_from = c(log2FC_combined, pvalue_combined, qvalue_combined),
              names_glue = "{Antibiotic}_{.value}")

# Save the combined dataframe to a new CSV file
write.csv(combined_wide, "combined_data_wide_format.csv", row.names = FALSE)

print("Combined data saved to 'combined_data_wide_format.csv'")
