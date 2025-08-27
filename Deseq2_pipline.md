# Bulk RNA-seq Analysis Workflow

## Step 1: Package Loading

### Convert `.tab` to `.csv` using Python

```python
import pandas as pd
import os

os.chdir("E:/GSE134031_DST120 (1).tab/")
data = pd.read_csv("GSE134031_DST120.tab", sep="\t")
data.to_csv("total_count_RNAseq.csv")
```

### Preprocessing with R (dplyr)

```R
library(DESeq2)
library(dplyr)

total_condition <- read.csv("condition.csv", header = TRUE, row.names = 1)
raw_count <- read.csv("total_count_RNAseq.csv", header = TRUE, row.names = 1)
raw_count2 <- raw_count[, -c(2,3)]

raw_count2_summarized <- raw_count2 %>%
  group_by(mgi_symbol) %>%
  summarize(across(where(is.numeric), ~sum(., na.rm = TRUE)))

raw_count2 <- as.data.frame(raw_count2_summarized)
row.names(raw_count2) <- raw_count2$mgi_symbol
raw_count2$mgi_symbol <- NULL
```

---

## Step 2: Exploratory Data Analysis





####condition
```powershell
PS E:\GSE134031_DST120 (1).tab> Get-Content condition.csv
,condition
IL1_M,TREM2_WT-young
IL10_M,TREM2_KO-young
IL11_M,TREM2_WT-old
IL12_M,TREM2_KO-old
IL13_M,TREM2_WT-young
IL14_M,TREM2_KO-young
IL15_M,TREM2_WT-old
IL16_M,TREM2_KO-old
IL17_M,TREM2_WT-young
IL18_M,TREM2_KO-young
IL19_M,TREM2_WT-old
IL2_M,TREM2_KO-young
IL20_M,TREM2_KO-old
IL21_M,TREM2_WT-young
IL22_M,TREM2_KO-young
IL23_M,TREM2_WT-old
IL24_M,TREM2_KO-old
IL25_M,TREM2_WT-young
IL26_M,TREM2_KO-young
IL27_M,TREM2_WT-old
IL28_M,TREM2_KO-old
IL3_M,TREM2_WT-old
IL4_M,TREM2_KO-old
IL5_M,TREM2_WT-young
IL6_M,TREM2_KO-young
IL7_M,TREM2_WT-old
IL8_M,TREM2_KO-old
IL9_M,TREM2_WT-young


```

### PCA and UMAP

```R
library(data.table)
library(ggplot2)
library(umap)



raw_count2 <- raw_count2[, row.names(total_condition)]
raw_count2 <- raw_count2[rowSums(raw_count2) > 0, ]

samples <- colnames(raw_count2)
condition_vector <- total_condition[samples, , drop = FALSE]$condition

log_counts <- log2(raw_count2 + 1)
log_counts_t <- t(log_counts)
log_counts_scaled <- scale(log_counts_t)

pca_result <- prcomp(log_counts_scaled, center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$Condition <- condition_vector

umap_result <- umap(log_counts_scaled)
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Condition <- condition_vector
```

---

## Step 3: Group-Specific Analysis

### Example: Old Mouse Group

```R
condition <- subset(total_condition, total_condition$condition %in% c("TREM2_WT-old", "TREM2_KO-old"))
raw_count2 <- raw_count2[, row.names(condition)]
raw_count2 <- raw_count2[rowSums(raw_count2) > 0, ]

samples <- colnames(raw_count2)
condition_vector <- condition[samples, , drop = FALSE]$condition

log_counts <- log2(raw_count2 + 1)
log_counts_t <- t(log_counts)
log_counts_scaled <- scale(log_counts_t)

pca_result <- prcomp(log_counts_scaled, center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$Condition <- condition_vector
```

### Example: Young Mouse Group

(Same workflow applied, using conditions `TREM2_WT-young` and `TREM2_KO-young`).

---

## Step 4: Differential Expression Analysis (DESeq2)

```R
library(DESeq2)

condition_df <- data.frame(condition = factor(condition$condition))
condition_df$condition <- relevel(condition_df$condition, ref = "TREM2_WT-old")
rownames(condition_df) <- rownames(condition)

raw_count_filtered <- raw_count2[, rownames(condition_df)]

dds <- DESeqDataSetFromMatrix(countData = raw_count_filtered,
                              colData = condition_df,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj), ]

write.csv(res, file = "DEG_TREM_KOvsWT.csv")
```

---

## Step 5: Visualization

### Volcano Plot

```R
library(EnhancedVolcano)

EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue',
    pCutoff = 0.05,
    FCcutoff = 1.5,
    title = 'TREM2_KO-old vs TREM2_WT-old',
    subtitle = 'DESeq2 results',
    pointSize = 1.0,
    labSize = 3.0)
```

### Non-Negative Matrix Factorization (NMF)

```R
library(NMF)

vst_data <- vst(dds, blind=FALSE)
vst_mat <- assay(vst_data)
vst_mat_top <- vst_mat[order(rowVars(vst_mat), decreasing = TRUE)[1:1000], ]

nmf_res <- nmf(vst_mat_top, rank=2:4, nrun=30, seed=123456)
plot(nmf_res)
```

---

## Step 6: Functional Enrichment Analysis

```R
library(clusterProfiler)
library(org.Mm.eg.db)

DEG_UP <- subset(res, res$log2FoldChange >= 1 & res$pvalue <= 0.05)
DEG_UP$GeneID <- row.names(DEG_UP)

filtered_DEG <- DEG_UP %>%
  dplyr::filter(!grepl("^ens", GeneID, ignore.case = TRUE) &
                !grepl("rik$", GeneID, ignore.case = TRUE) &
                !grepl("^rps", GeneID, ignore.case = TRUE) &
                !grepl("^rpl", GeneID, ignore.case = TRUE) &
                !grepl("^gm", GeneID, ignore.case = TRUE))

gene <- bitr(filtered_DEG$GeneID,
             fromType = "SYMBOL",
             toType = "ENTREZID",
             OrgDb = "org.Mm.eg.db")

ego_ALL <- enrichGO(gene = gene$ENTREZID,
                    OrgDb = org.Mm.eg.db,
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 1,
                    qvalueCutoff = 1,
                    readable = TRUE)
```

---

## Notes

* Outlier samples were removed iteratively (e.g., IL8\_M, IL12\_M, IL23\_M, IL27\_M, etc.) to improve PCA clustering.
* Both young and old groups were analyzed separately, followed by intersection analysis.
* Downstream analyses included GO and KEGG enrichment using **clusterProfiler**.
