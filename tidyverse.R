install.packages("tidyverse")
install.packages(
  c("arrow", "babynames", "curl", "duckdb", "gapminder", 
    "ggrepel", "ggridges", "ggthemes", "hexbin", "janitor", "Lahman", 
    "leaflet", "maps", "nycflights13", "openxlsx", "palmerpenguins", 
    "repurrrsive", "tidymodels", "writexl")
)

library(palmerpenguins)
library(ggthemes)
library(tidyverse)

# 读取 fpkm_tracking 文件
fpkm_data <- read_tsv("~/Library/CloudStorage/OneDrive-SaintLouisUniversity/BIOL_5930_shili/GSE73522_RAW/GSM1897222_CIN.genes.fpkm_tracking")

# 检查数据
glimpse(fpkm_data)  # 类似 str()
head(fpkm_data)     # 查看前几行
colnames(fpkm_data) # 检查列名

# 过滤表达量
filtered_data <- fpkm_data %>% 
  filter(FPKM > 1)  # 过滤FPKM大于1的基因

top_genes <- fpkm_data %>% 
  arrange(desc(FPKM)) %>% 
  slice_head(n = 10) # 或筛选top10高表达基因

# 计算统计量
summary_stats <- fpkm_data %>% 
  summarise(
    mean_FPKM = mean(FPKM, na.rm = TRUE),
    median_FPKM = median(FPKM, na.rm = TRUE),
    max_FPKM = max(FPKM, na.rm = TRUE)
  )

# 可视化
# fpkm分布直方图
ggplot(fpkm_data, aes(x = FPKM)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  scale_x_log10() +  # 如果FPKM跨度大，建议对数变换
  labs(title = "FPKM Distribution", x = "FPKM", y = "Count") +
  theme_minimal()

# 高表达基因柱状图
ggplot(top_genes, aes(x = reorder(tracking_id, -FPKM), y = FPKM)) +
  geom_col(fill = "tomato") +
  coord_flip() +  # 横向排列
  labs(title = "Top 10 Expressed Genes", x = "Gene", y = "FPKM") +
  theme_minimal()

# 计算不同基因的表达水平分布
fpkm_data %>%
  group_by(tracking_id) %>%
  summarise(mean_FPKM = mean(FPKM, na.rm = TRUE)) %>%
  arrange(desc(mean_FPKM)) %>%
  slice_head(n = 10)  # 获取平均表达量最高的前10个基因

ggplot(data = fpkm_data,mapping = aes())
