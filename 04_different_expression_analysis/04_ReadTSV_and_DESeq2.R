
# install DESeq2
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#  BiocManager::install("DESeq2")

library(DESeq2)
library(dplyr)
library(readr)   
library(tibble)
library(tidyverse)

# Set file path 
count_dir <- "~/Library/CloudStorage/OneDrive-SaintLouisUniversity/BIOL_5930_shili/abundance"  # 修改为你的实际路径
files <- list.files(count_dir, pattern = "_abundance.tsv$", full.names = TRUE)

# set sample name as colum name 
sample_names <- gsub("_abundance.tsv", "", basename(files))

# Initialization of the Merge Matrix 
# read file
df <- read.table(files[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
counts_matrix <- df %>% select(target_id, est_counts) %>%
  rename(!!sample_names[1] := est_counts)

# Cyclically read the remaining file and concatenate 
for (i in 2:length(files)) {
  temp <- read.table(files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  counts_matrix <- counts_matrix %>%
    left_join(temp %>% select(target_id, est_counts) %>%
                rename(!!sample_names[i] := est_counts),
              by = "target_id")
}

# set target_id as rownames
counts_matrix_final <- as.data.frame(counts_matrix)
rownames(counts_matrix_final) <- counts_matrix_final$target_id
counts_matrix_final <- counts_matrix_final[, -1]  # 去掉 target_id 列

# check results
head(counts_matrix_final)


###Save results
counts_matrix_final %>%
  rownames_to_column(var = "target_id") %>% 
  mutate(across(-target_id, round)) %>%        
  write_tsv("counts_matrix.tsv")

###Read tsv file

df <- read_tsv("counts_matrix.tsv") %>%
  column_to_rownames("target_id")
head(df)



###DESeq2 to analysis different gene expression
  

# read "counts_matrix.tsv"
counts_matrix <- read.table(
  file      = "counts_matrix.tsv",
  header    = TRUE,     # Have a header
  row.names = 1,        # Column 1 is the row name
  sep       = "\t"      # Separated by tabs
)

# check
head(counts_matrix)
dim(counts_matrix)  # Examine the number of rows and columns.

##read group infomation
sample_info <- read.table(
  file               = "abundance/sample.txt", 
  header             = TRUE,       # The first row is the column name
  sep                = "\t",       # tab-separated
  stringsAsFactors   = FALSE       # Don't convert strings to factors automatically
)
sample_info
str(sample_info)

#order
sample_info <- sample_info %>%
  arrange(match(sample, colnames(counts_matrix)))

# save file
#write_tsv(sample_info,file = "abundance/sample.txt")

#Ensure that the sample information is in consistent order with the columns of the count matrix.
sample_info <- sample_info[match(colnames(counts_matrix), sample_info$sample),]

###################################
# Create a DESeqDataSet object 
dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = sample_info,
                              design = ~ condition)
# Conducting DESeq2 Analysis
dds <- DESeq(dds)

# Obtaining Standardized Counts
normalized_counts <- counts(dds, normalized=TRUE)

# 1. save PCA picture 
vsd <- vst(dds, blind=TRUE)
pca_plot <- plotPCA(vsd, intgroup="condition") +
  theme_minimal() +
  ggtitle("PCA of Normalized Counts")
pca_plot

#ggsave("pca_plot.png", pca_plot, width=8, height=6, dpi=300)

# 2. Draw a heat map
library(pheatmap)

# Get the normalized count
normalized_counts <- counts(dds, normalized=TRUE)

# Select highly mutated genes (such as the top 10000 most mutated genes
var_genes <- order(rowVars(normalized_counts), decreasing=TRUE)[1:10000]
heatmap_data <- normalized_counts[var_genes,]
#heatmap_data <- normalized_counts

# Create sample comments and sort by condition
annotation_col <- data.frame(
  Condition = sample_info$condition,
  row.names = sample_info$sample
)


# order 
desired_order <- c("unfertilised_egg", "fertilised_egg", "16_cell", 
                   "initial_gastrula", "late_neurula", "mid_tailbud_II",
                   "late_tailbud_II", "larva")

# Sort sample information according to the specified order
sample_info_sorted <- sample_info %>%
  mutate(condition = factor(condition, levels = desired_order)) %>%
  arrange(condition)

# Rearrange the columns of heatmap_data according to the sorted samples

heatmap_data_sorted <- heatmap_data[, match(sample_info_sorted$sample, colnames(heatmap_data))]

# Create sample comments, make sure Condition is a factor and use the specified order
annotation_col <- data.frame(
  Condition = factor(sample_info$condition, levels = desired_order),
  row.names = sample_info$sample
)

# set color
#my_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# draw heatmap
pheatmap(heatmap_data,
         scale = "row",              
         clustering_distance_rows = "euclidean",
         cluster_cols = F,
         clustering_method = "complete",
         annotation_col = annotation_col,
         color = my_palette,
         show_rownames = FALSE,      
        #main = "Heatmap of Top 50 Variable Genes",
         filename = "heatmap.png",   
         width = 10,
         height = 8)

# Drawing heat map - Not correct column clustering 
pheatmap(heatmap_data_sorted,
         scale = "row",              
         cluster_cols = FALSE,       
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         annotation_col = annotation_col,
         color = my_palette,
         show_rownames = FALSE,
         #main = "Heatmap of Top 50 Variable Genes (Sorted by Condition)",
         filename = "heatmap_sorted.png",
         width = 10,
         height = 8)

################volcanoplot
#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

# Extract the differential expression results between the two conditions
#    The order of the contrast parameters is: variable name, condition for the numerator (baseline/control),
#    and condition for the denominator (treatment).
#    Note: The log2FoldChange in the results represents the log-transformed fold change of the baseline (first condition)
#    relative to the treatment (second condition).
res <- results(dds, contrast = c("condition", "initial_gastrula", "late_neurula"))


# Use the EnhancedVolcano package to plot the volcano plot
volcano_plot <- EnhancedVolcano(res,
                                lab = rownames(res),           # Gene labels
                                x = "log2FoldChange",           # Horizontal axis: log2 fold change
                                y = "pvalue",                   # Vertical axis: p-value (you can also use "padj" for adjusted p-value)
                                title = "initial_gastrula vs late_neurula",  # Plot title
                                pCutoff = 0.05,                 # p-value threshold
                                FCcutoff = 1.0,                 # log2 fold change threshold (adjust as needed)
                                pointSize = 2.5,                # Size of the points
                                labSize = 2.0)
# Display the volcano plot
volcano_plot

# Save the volcano plot as a PNG file with custom dimensions
ggsave("volcano_plot.png", plot = volcano_plot, width = 10, height = 16, dpi = 300)



# set color
keyvals <- ifelse(
  res$log2FoldChange > 1 & res$pvalue < 0.05, '#E3C6E0',
  ifelse(res$log2FoldChange < -1 & res$pvalue < 0.05, '#DBEDC5', '#BBBBBB'))

keyvals[is.na(keyvals)] <- 'black'


names(keyvals)[keyvals == '#E3C6E0'] <- 'up'
names(keyvals)[keyvals == '#DBEDC5'] <- 'down'
names(keyvals)[keyvals == '#BBBBBB'] <- 'non_sig'

# EnhancedVolcano
volcano_plot2 <- EnhancedVolcano(res,
               lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'initial_gastrula vs late_neurula',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 2,
                colCustom = keyvals,
                colAlpha = 0.75,
                legendLabSize = 15,
                legendIconSize = 5.0,
                arrowheads = FALSE,
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                border = 'partial',
                borderWidth = 1.5,
                labSize = 2.0,
                borderColour = 'black')
volcano_plot2
ggsave("volcano_plot.png", plot = volcano_plot2, width = 10, height = 10, dpi = 300)
