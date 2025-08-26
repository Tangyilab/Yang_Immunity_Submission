火山图作图 

```
library(org.Mm.eg.db)
library(tidyverse)

# 你要用的GO ID和对应注释
# 你要用的GO ID和对应注释
go_ids <- c("GO:0006935", "GO:0001935", "GO:0010573")
go_terms <- c("chemotaxis", "endothelial cell proliferation", "VEGF production")

# 自定义函数：通过 GO ID 获取 SYMBOL
get_symbols_by_go <- function(go_id) {
  # 获取ENTREZID
  entrez_ids <- get(go_id, org.Mm.egGO2ALLEGS)
  # 映射到SYMBOL
  symbols <- unlist(mget(entrez_ids, org.Mm.egSYMBOL, ifnotfound = NA))
  symbols <- symbols[!is.na(symbols)] # 去除NA
  return(symbols)
}

# 构建 data.frame：GeneID + GO注释名称
go_annot_df <- map2_dfr(go_ids, go_terms, function(go, term) {
  syms <- get_symbols_by_go(go)
  tibble(GeneID = syms, status = paste0("GO: ", term))
})


DEG_UP_young <- as.data.frame(res)
DEG_UP_young$GeneID <- rownames(DEG_UP_young)

# 筛掉非标准gene symbol的行
DEG_UP_young <- DEG_UP_young %>%
  filter(!grepl("^ens|rik$|^rps|^rpl|^rp23|^rp24|^mt-|^gm", GeneID, ignore.case = TRUE))

DEG_UP_young$status<-"NotSig"
DEG_UP_young <- DEG_UP_young %>%
  mutate(status_base = case_when(
    log2FoldChange > 0 & pvalue <= 0.05 ~ "Up",
    log2FoldChange < 0 & pvalue <= 0.05 ~ "Down",
    TRUE ~ "NotSig"
  )) %>%
  left_join(go_annot_df, by = "GeneID") %>%
  mutate(status = case_when(
    status_base %in% c("Up", "Down") & !is.na(status.y) ~ status.y,
    TRUE ~ status_base
  )) %>%
  dplyr::select(-status_base, -status.y)



library(ggplot2)
library(ggrepel)

# 去除NA
df <- na.omit(DEG_UP_young)

# 提取 GO 相关的基因用于标注
label_data <- df %>%
  filter(abs(log2FoldChange) >= 1 & 
         pvalue <= 0.05) %>%
  mutate(score = -log10(pvalue) * abs(log2FoldChange)) %>%
  arrange(desc(score)) %>%
  slice(1:90)
  
label_data  <- subset(label_data , status != "Up")  
label_data <- label_data %>%
  distinct(GeneID, .keep_all = TRUE)

df$status <- factor(df$status, levels = c("Up", "Down", "GO: chemotaxis", 
                                          "GO: endothelial cell proliferation", 
                                          "GO: VEGF production", "NotSig"))



# 自定义颜色
status_colors <- c(
  "Up" = "#E41A1C",
  "Down" = "#377EB8",
  "NotSig" = "grey70",
  "GO: chemotaxis" = "#4DAF4A",
  "GO: endothelial cell proliferation" = "#984EA3",
  "GO: VEGF production" = "#FF7F00"
)
df$status <- factor(df$status, levels = c("Up", "Down", "GO: chemotaxis", 
                                          "GO: endothelial cell proliferation", 
                                          "GO: VEGF production", "NotSig"))


 


highlight_genes <- c("Anxa1", "Ccl2", "Ccl24", "Ccl7", "Ccl8", "Ccr1", "Ccr2",
                     "Cxcl1", "Cxcr3", "Epha2", "Fpr2", "Gas6", "Hbegf", "Nr4a1",
                     "Il6", "Ptgs2", "Pdcl3", "Prkx", "Trem2")

# Step 2: 生成图
pdf("Valcano2.pdf", width = 8, height = 8)
ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue), color = status)) +
  geom_point(alpha = 0.7, size = 1.6) +
  scale_color_manual(values = status_colors) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot with Gene Labeling",
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Status"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "right") +
  scale_x_continuous(limits = c(-15, 15)) +
  scale_y_break(breaks = c(11.5, 33.5), space = 0.2, scales = 0.1, expand = c(0, 0)) +
  geom_text_repel(
    data = df %>% filter(GeneID %in% highlight_genes)%>% distinct(GeneID, .keep_all = TRUE),
    aes(label = GeneID),
    size = 4,
    max.overlaps = 100,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.size = 0.2
  )

dev.off()


# Step 2: 生成图
pdf("Valcano3.pdf", width = 8, height = 8)
ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue), color = status)) +
  geom_point(alpha = 0.7, size = 1.6) +
  scale_color_manual(values = status_colors) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot with Gene Labeling",
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Status"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "right") +
  scale_x_continuous(limits = c(-20, 20)) +
  scale_y_break(breaks = c(11.5, 33.5), space = 0.2, scales = 0.1, expand = c(0, 0)) +
  geom_text_repel(
    data = df %>% filter(GeneID %in% highlight_genes)%>% distinct(GeneID, .keep_all = TRUE),
    aes(label = GeneID),
    size = 4,
    max.overlaps = 100,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.size = 0.2
  )

dev.off()


```



<img src="C:\Users\dengy\AppData\Roaming\Typora\typora-user-images\image-20250415201614990.png" alt="image-20250415201614990" style="zoom:50%;" />

还凑合

制作一下基因热图 

```
total_condition <-read.csv("condition.csv",header = T,row.names = 1)
#读取总条件
raw_count <-read.csv("test.csv",header = T,row.names = 1)
#读取总矩阵
raw_count2<-raw_count[,-c(2,3)]
#删除无关信息
raw_count2_summarized <- raw_count2 %>%
  group_by(mgi_symbol) %>%
  summarize(across(where(is.numeric), ~sum(., na.rm = TRUE)))

# 将 tibble 转换为普通数据框
raw_count2 <- as.data.frame(raw_count2_summarized)
# 设置行名
row.names(raw_count2) <- raw_count2$mgi_symbol
# 移除原数据框中的 mgi_symbol 列（因为现在它是行名）
raw_count2$mgi_symbol <- NULL

total_condition <- total_condition %>%
  filter(!rownames(.) %in% c("IL8_M","IL3_M","IL20_M", "IL12_M", "IL23_M","IL27_M","IL6_M","IL2_M","IL17_M","IL13_M","IL1_M","IL10_M"))

raw_count2<-raw_count2[,row.names(total_condition)]
#过滤掉全是0的行
raw_count2 <- raw_count2[rowSums(raw_count2) > 0, ]




DEG_UP_young2 <- subset(DEG_UP_young,  log2FoldChange >=1 & pvalue <=0.05)
DEG_UP_young2 <- subset(DEG_UP_young2, status != "Up")

heatmap_gene <-unique(DEG_UP_young2$GeneID)
raw_count3<-raw_count2[heatmap_gene,]


library(pheatmap)

# 1. 转换数据为 matrix
mat <- as.matrix(raw_count3)

# 2. 获取列对应的 biological condition（从 total_condition）
col_annotation <- total_condition[match(colnames(mat), rownames(total_condition)), , drop = FALSE]

# 3. 设置排序顺序
condition_order <- c("TREM2_WT-young", "TREM2_KO-young", "TREM2_WT-old", "TREM2_KO-old")
col_annotation$condition <- factor(col_annotation$condition, levels = condition_order)

# 4. 按照condition排序列
mat <- mat[, order(col_annotation$condition)]

# 5. 重新匹配 annotation 排序
col_annotation <- col_annotation[order(col_annotation$condition), , drop = FALSE]

# 6. 作图
pdf("heatmap.pdf",width =7,height=8)
pheatmap(mat,
         annotation_col = col_annotation,
         cluster_cols = FALSE,  # 不聚类列
         cluster_rows = FALSE,   # 可以聚类基因
         show_colnames = FALSE, # 不显示列名（相当于去掉了x轴的标签）
         show_rownames = TRUE,  # 显示行名
         fontsize_row = 12,     # 增大y轴（行名）的字号到12
         fontsize_col = 8,      # 如果你仍需要显示列名，可以调整这个值；这里我们不显示列名
         scale = "row")      
dev.off()
```





