
library(dplyr)
library(tibble)
library(tidyverse)

# Set file path 设置文件路径
count_dir <- "~/Library/CloudStorage/OneDrive-SaintLouisUniversity/BIOL_5930_shili/abundance"  # 修改为你的实际路径
files <- list.files(count_dir, pattern = "_abundance.tsv$", full.names = TRUE)

# set sample name as colum name 提取样本名作为列名（从文件名中去掉路径和扩展名）
sample_names <- gsub("_abundance.tsv", "", basename(files))

# Initialization of the Merge Matrix 初始化合并矩阵
# read file
df <- read.table(files[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
counts_matrix <- df %>% select(target_id, est_counts) %>%
  rename(!!sample_names[1] := est_counts)

# Cyclically read the remaining file and concatenate 循环读取剩余文件并合并
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

library(tibble)
library(readr)

###Save results
counts_matrix_final %>%
  rownames_to_column(var = "target_id") %>% 
  mutate(across(-target_id, round)) %>%        # 对除 target_id 以外的所有列进行四舍五入
  write_tsv("counts_matrix.tsv")

###Read tsv file
library(readr)
library(dplyr)
library(tibble)

df <- read_tsv("counts_matrix.tsv") %>%
  column_to_rownames("target_id")
head(df)



###DESeq2

# install DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

library(DESeq2)
library(dplyr)
library(readr)   # If you wish to use functions such as `read_tsv`如果想用 read_tsv 等函数
library(tibble)  # If one wishes to use functions such as `rownames_to_column`, etc.如果想用 rownames_to_column 等函数

# read "counts_matrix.tsv"
counts_matrix <- read.table(
  file      = "counts_matrix.tsv",
  header    = TRUE,     # Have a header有表头
  row.names = 1,        # Column 1 is the row name (gene ID)第1列是行名(基因ID)
  sep       = "\t"      # Separated by tabs以制表符分隔
)

# check
head(counts_matrix)
dim(counts_matrix)  # Examine the number of rows and columns.看看行列数

##read group infomation
sample_info <- read.table(
  file               = "abundance/sample.txt", 
  header             = TRUE,       # The first row is the column name第一行是列名
  sep                = "\t",       # tab-separated以制表符分隔
  stringsAsFactors   = FALSE       # Don't convert strings to factors automatically不要把字符串自动转为因子
)
sample_info
str(sample_info)

#order
sample_info <- sample_info %>%
  arrange(match(sample, colnames(counts_matrix)))

# save file
#write_tsv(sample_info,file = "abundance/sample.txt")

#Ensure that the sample information is in consistent order with the columns of the count matrix.
# 确保样本信息与计数矩阵的列顺序一致
sample_info <- sample_info[match(colnames(counts_matrix), sample_info$sample),]

###################################
# Create a DESeqDataSet object 创建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = sample_info,
                              design = ~ condition)
# Conducting DESeq2 Analysis运行 DESeq2 分析
dds <- DESeq(dds)

# Obtaining Standardized Counts获取标准化后的计数
normalized_counts <- counts(dds, normalized=TRUE)

# 1. save PCA picture 保存 PCA 图
vsd <- vst(dds, blind=TRUE)
pca_plot <- plotPCA(vsd, intgroup="condition") +
  theme_minimal() +
  ggtitle("PCA of Normalized Counts")
pca_plot

# 保存为 PNG 文件
ggsave("pca_plot.png", pca_plot, width=8, height=6, dpi=300)

# 2. Draw a heat map绘制热图
library(pheatmap)

# Get the normalized count获取标准化计数
normalized_counts <- counts(dds, normalized=TRUE)

# Select highly mutated genes (such as the top 10000 most mutated genes选择高变异基因（例如前 10000 个最变异的基因）
var_genes <- order(rowVars(normalized_counts), decreasing=TRUE)[1:10000]
heatmap_data <- normalized_counts[var_genes,]
#heatmap_data <- normalized_counts

# Create sample comments and sort by condition创建样本注释并按 condition 排序
annotation_col <- data.frame(
  Condition = sample_info$condition,
  row.names = sample_info$sample
)


# order 定义您指定的顺序
desired_order <- c("unfertilised_egg", "fertilised_egg", "16_cell", 
                   "initial_gastrula", "late_neurula", "mid_tailbud_II",
                   "late_tailbud_II", "larva")

# Sort sample information according to the specified order按指定顺序对样本信息排序
sample_info_sorted <- sample_info %>%
  mutate(condition = factor(condition, levels = desired_order)) %>%
  arrange(condition)

# Rearrange the columns of heatmap_data according to the sorted samples
#根据排序后的样本顺序重新排列 heatmap_data 的列
heatmap_data_sorted <- heatmap_data[, match(sample_info_sorted$sample, colnames(heatmap_data))]

# Create sample comments, make sure Condition is a factor and use the specified order
#创建样本注释，确保 Condition 是因子并使用指定顺序
annotation_col <- data.frame(
  Condition = factor(sample_info$condition, levels = desired_order),
  row.names = sample_info$sample
)

# set color
my_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# draw heatmap
pheatmap(heatmap_data,
         scale = "row",              # 按行标准化
         clustering_distance_rows = "euclidean",
         cluster_cols = F,
         clustering_method = "complete",
         annotation_col = annotation_col,
         color = my_palette,
         show_rownames = FALSE,      # 如果基因名太多可以隐藏
        #main = "Heatmap of Top 50 Variable Genes",
         filename = "heatmap.png",   # 直接保存到文件
         width = 10,
         height = 8)

# Drawing heat map - Not correct column clustering 绘制热图 - 不对列聚类
pheatmap(heatmap_data_sorted,
         scale = "row",              # 按行标准化
         cluster_cols = FALSE,       # 不对列聚类
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         annotation_col = annotation_col,
         color = my_palette,
         show_rownames = FALSE,
         #main = "Heatmap of Top 50 Variable Genes (Sorted by Condition)",
         filename = "heatmap_sorted.png",
         width = 10,
         height = 8)
