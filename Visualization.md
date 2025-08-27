# Bulk RNA-seq Visualization Workflow

## Step 1: Volcano Plot


```R
library(org.Mm.eg.db)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(dplyr)

# Define GO IDs and corresponding annotation terms
go_ids <- c("GO:0006935", "GO:0001935", "GO:0010573")
go_terms <- c("chemotaxis", "endothelial cell proliferation", "VEGF production")

# Custom function: retrieve gene symbols from GO ID
get_symbols_by_go <- function(go_id) {
  entrez_ids <- get(go_id, org.Mm.egGO2ALLEGS)         # Get ENTREZ IDs
  symbols <- unlist(mget(entrez_ids, org.Mm.egSYMBOL, ifnotfound = NA)) # Map to SYMBOL
  symbols <- symbols[!is.na(symbols)]                  # Remove NA values
  return(symbols)
}

# Build annotation dataframe: GeneID + GO term
go_annot_df <- map2_dfr(go_ids, go_terms, function(go, term) {
  syms <- get_symbols_by_go(go)
  tibble(GeneID = syms, status = paste0("GO: ", term))
})

# Load DE results
DEG_UP_young <- as.data.frame(res)
DEG_UP_young$GeneID <- rownames(DEG_UP_young)

# Filter out non-standard gene symbols
DEG_UP_young <- DEG_UP_young %>%
  filter(!grepl("^ens|rik$|^rps|^rpl|^rp23|^rp24|^mt-|^gm", GeneID, ignore.case = TRUE))

# Assign status: Up, Down, NotSig
DEG_UP_young$status <- "NotSig"
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

# Remove NA rows
df <- na.omit(DEG_UP_young)

# Select significant genes for labeling
label_data <- df %>%
  filter(abs(log2FoldChange) >= 1 & pvalue <= 0.05) %>%
  mutate(score = -log10(pvalue) * abs(log2FoldChange)) %>%
  arrange(desc(score)) %>%
  slice(1:90)

label_data <- subset(label_data, status != "Up")  
label_data <- label_data %>%
  distinct(GeneID, .keep_all = TRUE)

# Define status factor order
df$status <- factor(df$status, levels = c(
  "Up", "Down", "GO: chemotaxis", 
  "GO: endothelial cell proliferation", 
  "GO: VEGF production", "NotSig"
))

# Custom color palette
status_colors <- c(
  "Up" = "#E41A1C",
  "Down" = "#377EB8",
  "NotSig" = "grey70",
  "GO: chemotaxis" = "#4DAF4A",
  "GO: endothelial cell proliferation" = "#984EA3",
  "GO: VEGF production" = "#FF7F00"
)

# Highlighted genes of interest
highlight_genes <- c("Anxa1", "Ccl2", "Ccl24", "Ccl7", "Ccl8", "Ccr1", "Ccr2",
                     "Cxcl1", "Cxcr3", "Epha2", "Fpr2", "Gas6", "Hbegf", "Nr4a1",
                     "Il6", "Ptgs2", "Pdcl3", "Prkx", "Trem2")

# Volcano Plot (version 1)
pdf("Volcano2.pdf", width = 8, height = 8)
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
    data = df %>% filter(GeneID %in% highlight_genes) %>% distinct(GeneID, .keep_all = TRUE),
    aes(label = GeneID),
    size = 4,
    max.overlaps = 100,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.size = 0.2
  )
dev.off()

# Volcano Plot (version 2, wider axis limits)
pdf("Volcano3.pdf", width = 8, height = 8)
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
    data = df %>% filter(GeneID %in% highlight_genes) %>% distinct(GeneID, .keep_all = TRUE),
    aes(label = GeneID),
    size = 4,
    max.overlaps = 100,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.size = 0.2
  )
dev.off()

```
# ===========================
## Step 2: Heatmap Visualization
# ===========================

```R
library(pheatmap)
library(tidyverse)

total_condition <- read.csv("condition.csv", header = TRUE, row.names = 1)
raw_count <- read.csv("total_count_RNAseq.csv", header = TRUE, row.names = 1)
raw_count2 <- raw_count[, -c(2,3)]

raw_count2_summarized <- raw_count2 %>%
  group_by(mgi_symbol) %>%
  summarize(across(where(is.numeric), ~sum(., na.rm = TRUE)))

raw_count2 <- as.data.frame(raw_count2_summarized)
row.names(raw_count2) <- raw_count2$mgi_symbol
raw_count2$mgi_symbol <- NULL

# Remove outlier samples
total_condition <- total_condition %>%
  filter(!rownames(.) %in% c("IL8_M","IL3_M","IL20_M", "IL12_M", "IL23_M","IL27_M",
                             "IL6_M","IL2_M","IL17_M","IL13_M","IL1_M","IL10_M"))

raw_count2 <- raw_count2[, row.names(total_condition)]
raw_count2 <- raw_count2[rowSums(raw_count2) > 0, ]

# Select significant DEGs for heatmap
DEG_UP_young2 <- subset(DEG_UP_young, log2FoldChange >= 1 & pvalue <= 0.05)
DEG_UP_young2 <- subset(DEG_UP_young2, status != "Up")

heatmap_gene <- unique(DEG_UP_young2$GeneID)
raw_count3 <- raw_count2[heatmap_gene, ]



# Prepare matrix and annotation
mat <- as.matrix(raw_count3)
col_annotation <- total_condition[match(colnames(mat), rownames(total_condition)), , drop = FALSE]

condition_order <- c("TREM2_WT-young", "TREM2_KO-young", "TREM2_WT-old", "TREM2_KO-old")
col_annotation$condition <- factor(col_annotation$condition, levels = condition_order)

mat <- mat[, order(col_annotation$condition)]
col_annotation <- col_annotation[order(col_annotation$condition), , drop = FALSE]

# Generate heatmap
pdf("heatmap.pdf", width = 7, height = 8)
pheatmap(mat,
         annotation_col = col_annotation,
         cluster_cols = FALSE,   # Do not cluster columns
         cluster_rows = FALSE,   # Optionally cluster rows (genes)
         show_colnames = FALSE,  # Hide column names
         show_rownames = TRUE,   # Show gene names
         fontsize_row = 12,
         fontsize_col = 8,
         scale = "row")
dev.off()
```
# ===========================
## Step 3: GSEA Visualization
# ===========================

```R

deg$external_gene_name = rownames(deg)
deg <- deg[deg$Status!='NotSig',]
DEG = deg[,c('external_gene_name','log2FoldChange','pvalue')]
head(DEG)
DEG <- DEG[order(DEG$log2FoldChange, decreasing = T),]
genelist <- DEG$log2FoldChange
names(genelist) <- DEG$external_gene_name

res_all <- GSEA(
  genelist,
  TERM2GENE = geneset,
  pvalueCutoff = 1,
  BPPARAM   = SerialParam()  # 串行执行
)

res_all_df = as.data.frame(res_all)
write.csv(as.data.frame(res_all),file = 'ALL_KO_WT_GSEA.csv')
deg$external_gene_name = rownames(deg)
deg <- deg[deg$Status!='NotSig',]
DEG = deg[,c('external_gene_name','log2FoldChange','pvalue')]
head(DEG)
DEG <- DEG[order(DEG$log2FoldChange, decreasing = T),]
genelist <- DEG$log2FoldChange
names(genelist) <- DEG$external_gene_name

res_all <- GSEA(
  genelist,
  TERM2GENE = geneset,
  pvalueCutoff = 1,
  BPPARAM   = SerialParam()  # 串行执行
)
res_all_df = as.data.frame(res_all)
write.csv(as.data.frame(res_all),file = 'ALL_KO_WT_GSEA.csv')

Terms = c(
    "GOBP_INTERLEUKIN_4_PRODUCTION",  #KO
    "GOBP_NEGATIVE_REGULATION_OF_VASCULATURE_DEVELOPMENT",  #KO
    "GOBP_POSITIVE_REGULATION_OF_CD4_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION",  #KO
    "GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS",  #KO
    "GOBP_POSITIVE_REGULATION_OF_REGULATORY_T_CELL_DIFFERENTIATION",   #WT
    "GOBP_NEGATIVE_REGULATION_OF_PRODUCTION_OF_MOLECULAR_MEDIATOR_OF_IMMUNE_RESPONSE",  #WT
    "GOBP_NEGATIVE_REGULATION_OF_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE",   #WT
    "GOBP_NEGATIVE_REGULATION_OF_PLATELET_ACTIVATION"
)
matched_rows <- res_all[grepl(paste(Terms, collapse = "|"), res_all$ID), ]
colnames(matched_rows)

devtools::install_github("junjunlab/GseaVis")
library(GseaVis)
head(res_all)
gseaNb(object = res_all,geneSetID ='GOBP_POSITIVE_REGULATION_OF_ACUTE_INFLAMMATORY_RESPONSE')

colors = colorRampPalette(c(
    "#376795", "#72BCD5", "#AADCE0", "white",
    "#FFE6B7", "#F7AA58", "#E76254"
  ))(100)
#08519C", "#A50F15"
genesets= c(
'GOBP_POSITIVE_REGULATION_OF_ACUTE_INFLAMMATORY_RESPONSE',
'GOBP_SPROUTING_ANGIOGENESIS',
'GOBP_POSITIVE_REGULATION_OF_IMMUNE_EFFECTOR_PROCESS'
,
'GOBP_NEGATIVE_REGULATION_OF_VASCULATURE_DEVELOPMENT'
)
p = gseaNb(object = res_all,geneSetID =genesets,subPlot = 2,
       #addPval = T,pvalX = 0.05,pvalY = 0.05,
       rmHt = F,htHeight = 0.5,htCol = c('#376795','#E76254') ,termWidth = 35,legend.position =c(0.8,0.7) ,
       curveCol = c("#84a494", "#dc6c4c", "#645cac","#d4c464"))
p +
  xlab("Trem2 KO - WT") +      # 修改 X 轴标签
  ylab("Enrichment Score") +    # 修改 Y 轴标签
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 如果需要旋转 X 轴标签，可以加上这个
p
ggsave("GSEA_Vis.pdf",width = 8,height = 5)

```

