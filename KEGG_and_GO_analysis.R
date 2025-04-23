
install.packages(c("tidyverse", "readxl", "janitor"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("clusterProfiler")

library(tidyverse)       
library(readxl)          
library(janitor)         # clean_names()
library(clusterProfiler) # enricher, barplot, dotplot

# 2. Read in the data and clean column names
fn <- "~/Library/CloudStorage/OneDrive-SaintLouisUniversity/BIOL_5930_shili/eggnog_mapper_results/out.emapper.annotations.xlsx"
emapper <- read_excel(fn, sheet = 1, skip = 2) %>% 
  clean_names()

# (Optional) Check the column names
print(colnames(emapper))

# 3. Rename key columns for clarity
emapper <- emapper %>%
  rename(
    gene   = query,
    go_raw = g_os,
    ko_raw = kegg_ko
  )


# 4. Split GO/KEGG annotations into long format
go_df <- emapper %>%
  filter(!is.na(go_raw) & go_raw != "-") %>%
  separate_rows(go_raw, sep = ",") %>%
  select(gene, GO = go_raw)

kegg_df <- emapper %>%
  filter(!is.na(ko_raw) & ko_raw != "-") %>%
  separate_rows(ko_raw, sep = ",") %>%
  select(gene, KO = ko_raw)

# 5. Define the background gene set
bg_genes <- unique(emapper$gene)

# 6. Build TERM2GENE mapping tables
go2gene   <- go_df   %>% rename(term = GO) %>% select(term, gene)
kegg2gene <- kegg_df %>% rename(term = KO) %>% select(term, gene)

# 7. Perform enrichment analysis
go_enrich <- enricher(
  gene         = unique(go_df$gene),
  TERM2GENE    = go2gene,
  universe     = bg_genes,
  pvalueCutoff = 0.05,
  pAdjustMethod= "BH"
)

kegg_enrich <- enricher(
  gene         = unique(kegg_df$gene),
  TERM2GENE    = kegg2gene,
  universe     = bg_genes,
  pvalueCutoff = 0.05,
  pAdjustMethod= "BH"
)

# 8. Visualization

# Assume go_enrich has been calculated above
# Extract GO enrichment results as a data.frame
res_go <- if (!is.null(go_enrich)) as.data.frame(go_enrich) else NULL

# Check and visualize GO enrichment or top frequent GO terms
if (is.null(res_go) || nrow(res_go) == 0) {
  message("⚠️ No significant GO enrichment (p<0.05); showing top 20 frequent GO terms")
  
  go_counts <- go_df %>%
    count(GO, name = "Count") %>%
    arrange(desc(Count)) %>%
    slice_head(n = 20)
  
  ggplot(go_counts, aes(x = Count, y = reorder(GO, Count))) +
    geom_col() +
    labs(title = "Top20 GO Term Frequency", x = "Gene Count", y = "GO Term") +
    theme_minimal()
  
} else {
  # Only plot barplot/dotplot if enrichment results exist
  barplot(go_enrich, showCategory = 20, title = "GO Term Enrichment")
  dotplot(go_enrich, showCategory = 20, title = "GO Term Dotplot")
}

# Plot the GO frequency barplot with a blue fill
ggplot(go_counts, aes(x = Count, y = reorder(GO, Count))) +
  geom_col(fill = "steelblue") +
  labs(title = "Top20 GO Term Frequency", x = "Gene Count", y = "GO Term") +
  theme_minimal()


### KEGG analysis
# Convert KEGG enrichment results to a data.frame
res_kegg <- if (!is.null(kegg_enrich)) as.data.frame(kegg_enrich) else NULL

# Check and visualize KEGG enrichment or top frequent KO terms
if (is.null(res_kegg) || nrow(res_kegg) == 0) {
  message(" #No significant KEGG enrichment (p<0.05); displaying top 20 frequent KO terms")
  
  kegg_counts <- kegg_df %>%
    count(KO, name = "Count") %>%
    arrange(desc(Count)) %>%
    slice_head(n = 20)
  
  ggplot(kegg_counts, aes(x = Count, y = reorder(KO, Count))) +
    geom_col() +
    labs(
      title = "Top20 KEGG KO Frequency",
      x = "Gene Count",
      y = "KO term"
    ) +
    theme_minimal()
  
} else {
  # If enrichment results exist, plot barplot and dotplot
  barplot(kegg_enrich, showCategory = 20, title = "KEGG Pathway Enrichment")
  dotplot(kegg_enrich, showCategory = 20, title = "KEGG Pathway Dotplot")
}

# Plot the KEGG frequency barplot with an orange fill
ggplot(kegg_counts, aes(x = Count, y = reorder(KO, Count))) +
  geom_col(fill = "tomato") +
  labs(title = "Top20 KEGG KO Frequency", x = "Gene Count", y = "KO term") +
  theme_minimal()

