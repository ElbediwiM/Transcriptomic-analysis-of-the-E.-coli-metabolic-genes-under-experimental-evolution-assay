library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(ggrepel)
library(circlize)

# Load the data
data <- read.csv("combined_data_wide_format.csv")

# View the first few rows
head(data)


# Reshape data for volcano plot
volcano_data <- data %>%
  pivot_longer(cols = starts_with("AMP") | starts_with("TET") | starts_with("KAN") | 
                 starts_with("CHL") | starts_with("CPR") | starts_with("DOX") | 
                 starts_with("ERY") | starts_with("FOX") | starts_with("TRM") | 
                 starts_with("NAL") | starts_with("NIT") | starts_with("TOB"),
               names_to = c("Antibiotic", "Metric"),
               names_sep = "_") %>%
  pivot_wider(names_from = Metric, values_from = value)

# Check the column names in the reshaped data
print(colnames(volcano_data))


# Create a column for labeling significant genes
volcano_data <- volcano_data %>%
  mutate(Significant = ifelse(abs(log2FC) > 3.5 & pvalue < 0.05, Genes, NA))
# Define manual colors for antibiotics
antibiotic_colors <- c(
  "AMP" = "red",       # Red for AMP
  "TET" = "blue",      # Blue for TET
  "KAN" = "green",     # Green for KAN
  "CHL" = "purple",    # Purple for CHL
  "CPR" = "orange",    # Orange for CPR
  "DOX" = "brown",     # Brown for DOX
  "ERY" = "pink",      # Pink for ERY
  "FOX" = "cyan",      # Cyan for FOX
  "TRM" = "magenta",   # Magenta for TRM
  "NAL" = "gray",      # Gray for NAL
  "NIT" = "darkgreen", # Dark green for NIT
  "TOB" = "black"       # Gold for TOB
)

# Plot volcano plot with labels
p <- ggplot(volcano_data, aes(x = log2FC, y = -log10(pvalue), color = Antibiotic)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(aes(label = Significant, colour = Antibiotic),  # Label significant genes
                  box.padding = 0.5,          # Adjust spacing
                  max.overlaps = Inf,         # Allow unlimited overlaps
                  size = 4,                   # Label text size
                  segment.color = "gray50") +  # Line color for labels
  theme_minimal() +
  labs(title = "",
       x = "Log2 Fold Change",
       y = "-Log10(p-value)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Significance threshold
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "blue") +    # Fold change thresholds
  scale_color_manual(values = antibiotic_colors) +
  theme(legend.position = "right",  # Position of the legend
        legend.text = element_text(size = 14, face = "bold"),  # Bold and enlarge legend text
        legend.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"))  # Bold and enlarge legend title


# Save the volcano plot as a high-resolution image
ggsave("volcano_plot_high_res.png", plot = p, 
       width = 12, height = 8, dpi = 600, units = "in")

# Save the reshaped data as a CSV file
write.csv(volcano_data, "volcano_plot_data.csv", row.names = FALSE)

print("Volcano plot saved as 'volcano_plot_high_res.png'")
print("Volcano plot data saved as 'volcano_plot_data.csv'")







#####################################

library(ComplexHeatmap)
library(circlize)  # Required for colorRamp2
library(dplyr)
library(tidyr)

# Load the data
data <- read.csv("combined_data_wide_format.csv")

# Reshape data to include log2FC and p-value
volcano_data <- data %>%
  pivot_longer(cols = starts_with("AMP") | starts_with("TET") | starts_with("KAN") | 
                 starts_with("CHL") | starts_with("CPR") | starts_with("DOX") | 
                 starts_with("ERY") | starts_with("FOX") | starts_with("TRM") | 
                 starts_with("NAL") | starts_with("NIT") | starts_with("TOB"),
               names_to = c("Antibiotic", "Metric"),
               names_sep = "_") %>%
  pivot_wider(names_from = Metric, values_from = value)

# Identify significant genes (|log2FC| > 2 and p-value < 0.05 for any antibiotic)
significant_genes <- volcano_data %>%
  group_by(Genes) %>%
  filter(any(abs(log2FC) > 2 & pvalue < 0.05)) %>%
  ungroup()

# Reshape back to wide format to include all antibiotics
significant_wide <- significant_genes %>%
  select(Genes, Antibiotic, log2FC) %>%
  pivot_wider(names_from = Antibiotic, values_from = log2FC)

# Replace NA with 0
significant_wide <- significant_wide %>%
  replace(is.na(.), 0)

# Convert to a matrix
log2fc_matrix <- as.matrix(significant_wide[, -1])  # Exclude the "Genes" column
rownames(log2fc_matrix) <- significant_wide$Genes   # Set gene names as row names

# Sort the matrix by antibiotic (columns)
log2fc_matrix <- log2fc_matrix[, order(colnames(log2fc_matrix))]

# Define the color scale using colorRamp2
col_fun <- colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))

# Define colors for antibiotics (same as in the volcano plot)
antibiotic_colors <- c(
  "AMP" = "red",       # Red for AMP
  "TET" = "blue",      # Blue for TET
  "KAN" = "green",     # Green for KAN
  "CHL" = "purple",    # Purple for CHL
  "CPR" = "orange",    # Orange for CPR
  "DOX" = "brown",     # Brown for DOX
  "ERY" = "pink",      # Pink for ERY
  "FOX" = "cyan",      # Cyan for FOX
  "TRM" = "magenta",   # Magenta for TRM
  "NAL" = "gray",      # Gray for NAL
  "NIT" = "darkgreen", # Dark green for NIT
  "TOB" = "black"       # Gold for TOB
)

# Create a column annotation (colored strips on top of columns)
column_ha <- HeatmapAnnotation(
  Antibiotic = colnames(log2fc_matrix),  # Antibiotic names
  col = list(Antibiotic = antibiotic_colors),  # Use the defined colors
  annotation_name_side = "left",  # Place annotation names on the left
  show_legend = TRUE  # Show legend for the annotation
)


# Create the heatmap
heatmap <- Heatmap(
  log2fc_matrix,
  name = "Log2FC",  # Name of the legend
  column_title = "",  # Title for columns
  row_title = "Genes",           # Title for rows
  row_names_side = "right",      # Place row names on the right
  column_names_rot = 45,         # Rotate column names for better readability
  cluster_columns = TRUE,       # Do not cluster columns (keep sorted by antibiotic)
  cluster_rows = TRUE,           # Cluster rows (genes)
  show_row_names = TRUE,         # Show row names (genes)
  show_column_names = F,      # Show column names (antibiotics)
  col = col_fun,                 # Use the defined color scale
  top_annotation = column_ha,    # Add column annotations
  heatmap_legend_param = list(
    title = "Log2FC",            # Legend title
    at = c(-5, 0, 5),            # Legend breaks
    labels = c("-5", "0", "5")   # Legend labels
  )
)

# Draw the heatmap
draw(heatmap)

# Save the heatmap as a high-resolution PNG file
png("heatmap_significant_genes_high_res.png", width = 10, height = 10, units = "in", res = 300)
draw(heatmap)
dev.off()

print("Heatmap saved as 'heatmap_significant_genes_high_res.png'")


##########################################################################################

