BulkRNAseq 分析流程

加载包

使用pyton将tab转为csv

```
import pandas as pd

import os

os.chdir("E:/GSE134031_DST120 (1).tab/")

data = pd.read_csv("GSE134031_DST120.tab",sep="\t")

data.to_csv("test.csv")
```



接下来首先使用dplyr预先处理数据



```

library(DESeq2)
library(dplyr)

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
```



先提取信息，跑PCA和UMAP看看数据分布



```
library(data.table)
library(ggplot2)
library(umap)


raw_count2<-raw_count2[,row.names(total_condition)]
#过滤掉全是0的行
raw_count2 <- raw_count2[rowSums(raw_count2) > 0, ]


samples <- colnames(raw_count2)
condition_vector <- total_condition[samples, , drop = FALSE]$condition

# log2(TPM+1) 风格转换，避免对0取对数
log_counts <- log2(raw_count2 + 1)

# 转置为样本 x 基因
log_counts_t <- t(log_counts)

# 进行PCA前先中心化标准化
log_counts_scaled <- scale(log_counts_t)

pca_result <- prcomp(log_counts_scaled, center = TRUE, scale. = TRUE)

# 提取前两主成分
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$Condition <- condition_vector

# 可视化
ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA of Samples")
  
set.seed(123)  # 保证可重复
umap_result <- umap(log_counts_scaled)

umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Condition <- condition_vector

# 可视化
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Condition)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "UMAP of Samples")  
  
  

```



![image-20250414161915077](C:\Users\dengy\AppData\Roaming\Typora\typora-user-images\image-20250414161915077.png)

<img src="C:\Users\dengy\AppData\Roaming\Typora\typora-user-images\image-20250414162011388.png" alt="image-20250414162011388" style="zoom:80%;" />





看来数据不能全部放在一起分析，需要切割成年轻和老年小鼠 



先看2mon小鼠的差异



```
library(data.table)
library(ggplot2)
library(umap)

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


condition <- subset(total_condition,total_condition$condition %in% c("TREM2_WT-young", "TREM2_KO-young"))

raw_count2<-raw_count2[,row.names(condition)]
#过滤掉全是0的行
raw_count2 <- raw_count2[rowSums(raw_count2) > 0, ]


samples <- colnames(raw_count2)
condition_vector <- condition[samples, , drop = FALSE]$condition

# log2(TPM+1) 风格转换，避免对0取对数
log_counts <- log2(raw_count2 + 1)

# 转置为样本 x 基因
log_counts_t <- t(log_counts)

# 进行PCA前先中心化标准化
log_counts_scaled <- scale(log_counts_t)

pca_result <- prcomp(log_counts_scaled, center = TRUE, scale. = TRUE)

# 提取前两主成分
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$Condition <- condition_vector

# 可视化
ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA of Samples")
```

![image-20250414162909326](C:\Users\dengy\AppData\Roaming\Typora\typora-user-images\image-20250414162909326.png)



```
pca_df$SampleName <- rownames(pca_df)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) + # 绘制点
  geom_text(aes(label = rownames(pca_df)), vjust = -0.5, hjust = 0.5) + # 添加文本标签
  theme_minimal() +
  labs(title = "PCA of Samples with Sample Names") +
  theme(legend.position = "right") # 根据需要调整主题
```

<img src="C:\Users\dengy\AppData\Roaming\Typora\typora-user-images\image-20250414222013375.png" alt="image-20250414222013375" style="zoom:50%;" />



很散

看看old组 



```
library(data.table)
library(ggplot2)
library(umap)

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


condition <- subset(total_condition,total_condition$condition %in% c("TREM2_WT-old", "TREM2_KO-old"))

raw_count2<-raw_count2[,row.names(condition)]
#过滤掉全是0的行
raw_count2 <- raw_count2[rowSums(raw_count2) > 0, ]


samples <- colnames(raw_count2)
condition_vector <- condition[samples, , drop = FALSE]$condition

# log2(TPM+1) 风格转换，避免对0取对数
log_counts <- log2(raw_count2 + 1)

# 转置为样本 x 基因
log_counts_t <- t(log_counts)

# 进行PCA前先中心化标准化
log_counts_scaled <- scale(log_counts_t)

pca_result <- prcomp(log_counts_scaled, center = TRUE, scale. = TRUE)

# 提取前两主成分
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$Condition <- condition_vector

# 可视化
ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA of Samples")
```

![image-20250414164304455](C:\Users\dengy\AppData\Roaming\Typora\typora-user-images\image-20250414164304455.png)



似乎也很散，加一个标签，看看样本信息



```
pca_df$SampleName <- rownames(pca_df)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) + # 绘制点
  geom_text(aes(label = rownames(pca_df)), vjust = -0.5, hjust = 0.5) + # 添加文本标签
  theme_minimal() +
  labs(title = "PCA of Samples with Sample Names") +
  theme(legend.position = "right") # 根据需要调整主题
```



![image-20250414164553450](C:\Users\dengy\AppData\Roaming\Typora\typora-user-images\image-20250414164553450.png)

各删2个离群值，变成 5对5

IL8_M IL12_M IL23_M IL27_M  IL3_M

```
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


condition <- subset(total_condition,total_condition$condition %in% c("TREM2_WT-old", "TREM2_KO-old"))
###挑选老年组

condition <- condition %>%
  filter(!rownames(.) %in% c("IL8_M","IL3_M","IL20_M", "IL12_M", "IL23_M", "IL27_M"))
#删除  IL8_M IL12_M IL23_M IL27_M 
raw_count2<-raw_count2[,row.names(condition)]
#过滤掉全是0的行
raw_count2 <- raw_count2[rowSums(raw_count2) > 0, ]


samples <- colnames(raw_count2)
condition_vector <- condition[samples, , drop = FALSE]$condition

# log2(TPM+1) 风格转换，避免对0取对数
log_counts <- log2(raw_count2 + 1)

# 转置为样本 x 基因
log_counts_t <- t(log_counts)

# 进行PCA前先中心化标准化
log_counts_scaled <- scale(log_counts_t)

pca_result <- prcomp(log_counts_scaled, center = TRUE, scale. = TRUE)

# 提取前两主成分
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$Condition <- condition_vector


# 可视化
ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA of Samples
  
  
# 加载必要的库
library(ggplot2)

# 创建PDF文件
pdf("PCA_old.pdf", width = 5.5, height = 4)

# 假设pca_df包含列PC1, PC2和Condition
ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) + # 设置点的大小
  theme_bw() + # 使用黑白主题
  labs(title = "PCA of young old") + # 添加标题
  # 手动设置颜色
  scale_color_manual(values = c("#D1352B", "#B383B9")) +
  # 调整坐标轴文本和其他元素
  theme(
    axis.title.x = element_text(color = "black", size = 14), # x轴标签字体颜色和大小
    axis.title.y = element_text(color = "black", size = 14), # y轴标签字体颜色和大小
    axis.text.x = element_text(color = "black", size = 12), # x轴文本字体颜色和大小
    axis.text.y = element_text(color = "black", size = 12), # y轴文本字体颜色和大小
    plot.title = element_text(hjust = 0.5) # 标题居中
  )

# 关闭设备以完成PDF输出
dev.off()  
  
  
```

![image-20250414175824253](C:\Users\dengy\AppData\Roaming\Typora\typora-user-images\image-20250414175824253.png)

加个标签的版本

```
pca_df$SampleName <- rownames(pca_df)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) + # 绘制点
  geom_text(aes(label = rownames(pca_df)), vjust = -0.5, hjust = 0.5) + # 添加文本标签
  theme_minimal() +
  labs(title = "PCA of Samples with Sample Names") +
  theme(legend.position = "right") # 根据需要调整主题
```

![image-20250414175806004](C:\Users\dengy\AppData\Roaming\Typora\typora-user-images\image-20250414175806004.png)

开始走deseq2



```
library(DESeq2)


condition_df <- data.frame(condition = factor(condition$condition))
condition_df$condition <- relevel(condition_df$condition, ref = "TREM2_WT-old")
rownames(condition_df) <- rownames(condition)

# 确保 count matrix 列顺序和 condition 行顺序一致
raw_count_filtered <- raw_count2[, rownames(condition_df)]

# 创建 DESeq2 对象
dds <- DESeqDataSetFromMatrix(countData = raw_count_filtered,
                              colData = condition_df,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)


resOrdered <- res[order(res$padj), ]
head(resOrdered)

log2 fold change (MLE): condition TREM2 WT.old vs TREM2 KO.old 
Wald test p-value: condition TREM2 WT.old vs TREM2 KO.old 
DataFrame with 6 rows and 6 columns
        baseMean log2FoldChange     lfcSE      stat      pvalue        padj
       <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
Spp1     90.2532        4.61065  0.339811  13.56829 6.17436e-42 5.02778e-38
Lpl     124.3031        2.94339  0.222206  13.24619 4.74672e-40 1.93263e-36
Trem2   156.4916        2.28839  0.182500  12.53908 4.56269e-36 1.23847e-32
Sep15   142.9031       -1.52816  0.222898  -6.85588 7.08752e-12 1.44284e-08
Clec7a   57.4405        1.62440  0.238535   6.80988 9.76795e-12 1.59081e-08
F13a1   145.3777       -1.13036  0.169778  -6.65787 2.77831e-11 3.77063e-08

write.csv(res,file="DEG_TREM_KOvsWT.csv")




```



做一下火山图



```
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



<img src="C:\Users\dengy\AppData\Roaming\Typora\typora-user-images\image-20250414203915147.png" alt="image-20250414203915147" style="zoom:50%;" />



感觉有点丑

后续再调，走一下NMF看看

```
library(DESeq2)
library(NMF)

# 获取 vst 转换表达矩阵
vst_data <- vst(dds, blind=FALSE)
vst_mat <- assay(vst_data)

# 可选：选 top variable genes
vst_mat_top <- vst_mat[order(rowVars(vst_mat), decreasing = TRUE)[1:1000], ]

# 运行 NMF，设定聚类数 rank = 2~4（你也可以探索更多）
nmf_res <- nmf(vst_mat_top, rank=2:4, nrun=30, seed=123456)

# 查看聚类质量（cophenetic系数）
plot(nmf_res)

# 选出最佳的聚类数（如 rank = 2）
nmf_res <- nmf(vst_mat_top, rank=2, nrun=30, seed=123456)

# 展示一致性矩阵图
consensusmap(nmf_res)
```



<img src="C:\Users\dengy\AppData\Roaming\Typora\typora-user-images\image-20250414210949945.png" alt="image-20250414210949945" style="zoom:50%;" />



感觉意义不大，就是一个聚类

做个差异富集看看

```

DEG_UP <-subset(res,res$log2FoldChange>=1 & res$pvalue <=0.05)
head(DEG_UP)
log2 fold change (MLE): condition TREM2 KO.old vs TREM2 WT.old 
Wald test p-value: condition TREM2 KO.old vs TREM2 WT.old 
DataFrame with 6 rows and 6 columns
        baseMean log2FoldChange     lfcSE      stat      pvalue        padj
       <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
Abracl   24.7946        1.44267  0.326777   4.41484 1.01084e-05 2.16613e-03
Anxa1    11.8257        2.74638  0.680380   4.03653 5.42469e-05 8.18023e-03
Ccr1     32.0591        1.24887  0.313723   3.98079 6.86861e-05 9.81247e-03
Ccr2     56.6624        1.13223  0.252888   4.47719 7.56313e-06 1.71074e-03
Cd209a   23.3102        2.54370  0.421850   6.02988 1.64086e-09 1.48461e-06
Cd209f   44.5710        1.06053  0.272000   3.89902 9.65835e-05 1.15659e-02

#去除非编码基因 
DEG_UP$GeneID <-row.names(DEG_UP)
DEG_UP <-as.data.frame(DEG_UP)
filtered_DEG <- DEG_UP %>%
  dplyr::filter(!grepl("^ens", GeneID, ignore.case = TRUE) & 
                !grepl("rik$", GeneID, ignore.case = TRUE) &
                !grepl("^rps", GeneID, ignore.case = TRUE) &
                !grepl("^rpl", GeneID, ignore.case = TRUE) &
                !grepl("^gm", GeneID, ignore.case = TRUE))

                
#开始走GO和KEGG分析 
library(clusterProfiler)
library(org.Mm.eg.db)

gene = bitr(filtered_DEG$GeneID, #数据集
fromType="SYMBOL", #输入为SYMBOL格式
toType="ENTREZID",  # 转为ENTERZID格式
OrgDb="org.Mm.eg.db") #小鼠数据库 

head(gene)
  SYMBOL ENTREZID
2 Abca13   268379
3 Abracl    73112
4 Adam19    11492
5  Adam9    11502
6  Anxa1    16952
7  Apba3    57267


ego_ALL <- enrichGO(gene = gene$ENTREZID,
                OrgDb = org.Mm.eg.db, #没有organism="human"，改为OrgDb=org.Hs.eg.db
				#keytype = 'ENSEMBL',
                ont = "ALL", #也可以是 CC  BP  MF中的一种
                pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                pvalueCutoff = 1, #P值会过滤掉很多，可以全部输出
                qvalueCutoff = 1,
				readable = TRUE) #Gene ID 转成gene Symbol ，易读






                
```

似乎 有效的信息很少 我们试试把年轻组和年老组取个交集看看

走一下原始的rawdata看看

```
library(DESeq2)
library(NMF)

# raw_count2: 行是基因，列是样本

# 过滤掉表达很低的基因（比如在多数样本中都是0）
keep_genes <- rowSums(raw_count2 > 1) >= 2
expr_mat <- raw_count2[keep_genes, ]

#log转换一下
log_expr <- log2(expr_mat + 1)

# 设定重复次数（越高越稳定，但越慢）
nmf_res <- nmf(log_expr, 2:10, nrun = 30, seed = 123)

# 可视化 rank survey 图
plot(nmf_res)

##rank3是最优解 
best_fit <- nmf(log_expr, rank = 3, nrun = 30, seed = 123)

consensusmap(best_fit)

```



<img src="C:\Users\dengy\AppData\Roaming\Typora\typora-user-images\image-20250415094704853.png" alt="image-20250415094704853" style="zoom:50%;" />

一样的散 看来年轻组不适合分析，尝试去除一下  IL6_M IL2_M IL17_M  IL13_M

```
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


condition <- subset(total_condition,total_condition$condition %in% c("TREM2_WT-young", "TREM2_KO-young"))
###挑选老年组

condition <- condition %>%
  filter(!rownames(.) %in% c("IL6_M","IL2_M","IL17_M","IL13_M"))
#删除  IL8_M IL12_M IL23_M IL27_M 
raw_count2<-raw_count2[,row.names(condition)]
#过滤掉全是0的行
raw_count2 <- raw_count2[rowSums(raw_count2) > 0, ]


samples <- colnames(raw_count2)
condition_vector <- condition[samples, , drop = FALSE]$condition

# log2(TPM+1) 风格转换，避免对0取对数
log_counts <- log2(raw_count2 + 1)

# 转置为样本 x 基因
log_counts_t <- t(log_counts)

# 进行PCA前先中心化标准化 
log_counts_scaled <- scale(log_counts_t)

pca_result <- prcomp(log_counts_scaled, center = TRUE, scale. = TRUE)

# 提取前两主成分
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$Condition <- condition_vector


# 可视化
ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA of Samples")
```

<img src="C:\Users\dengy\AppData\Roaming\Typora\typora-user-images\image-20250415095619534.png" alt="image-20250415095619534" style="zoom:50%;" />

```
pca_df$SampleName <- rownames(pca_df)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) + # 绘制点
  geom_text(aes(label = rownames(pca_df)), vjust = -0.5, hjust = 0.5) + # 添加文本标签
  theme_minimal() +
  labs(title = "PCA of Samples with Sample Names") +
  theme(legend.position = "right") # 根据需要调整主题
```

<img src="C:\Users\dengy\AppData\Roaming\Typora\typora-user-images\image-20250415095706465.png" alt="image-20250415095706465" style="zoom:50%;" />

再走一下NMF

```
library(NMF)

# raw_count2: 行是基因，列是样本

# 过滤掉表达很低的基因（比如在多数样本中都是0）
keep_genes <- rowSums(raw_count2 > 1) >= 2
expr_mat <- raw_count2[keep_genes, ]

#log转换一下
log_expr <- log2(expr_mat + 1)

# 设定重复次数（越高越稳定，但越慢）
nmf_res <- nmf(log_expr, 2:10, nrun = 30, seed = 123)

# 可视化 rank survey 图
plot(nmf_res)

##rank3是最优解 
best_fit <- nmf(log_expr, rank = 3, nrun = 30, seed = 123)

consensusmap(best_fit)
```

<img src="C:\Users\dengy\AppData\Roaming\Typora\typora-user-images\image-20250415100706189.png" alt="image-20250415100706189" style="zoom:50%;" />

再删掉IL1_M IL10_M 

```
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


condition <- subset(total_condition,total_condition$condition %in% c("TREM2_WT-young", "TREM2_KO-young"))
###挑选老年组

condition <- condition %>%
  filter(!rownames(.) %in% c("IL6_M","IL2_M","IL17_M","IL13_M","IL1_M","IL10_M"))
#删除  IL8_M IL12_M IL23_M IL27_M 
raw_count2<-raw_count2[,row.names(condition)]
#过滤掉全是0的行
raw_count2 <- raw_count2[rowSums(raw_count2) > 0, ]


samples <- colnames(raw_count2)
condition_vector <- condition[samples, , drop = FALSE]$condition

# log2(TPM+1) 风格转换，避免对0取对数
log_counts <- log2(raw_count2 + 1)

# 转置为样本 x 基因
log_counts_t <- t(log_counts)

# 进行PCA前先中心化标准化
log_counts_scaled <- scale(log_counts_t)

pca_result <- prcomp(log_counts_scaled, center = TRUE, scale. = TRUE)

# 提取前两主成分
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$Condition <- condition_vector


# 可视化
pdf("PCA_young.pdf",width =6,height =4)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) +
  theme_bw() +
  labs(title = "PCA of young group")
dev.off()  
```

<img src="C:\Users\dengy\AppData\Roaming\Typora\typora-user-images\image-20250415100845668.png" alt="image-20250415100845668" style="zoom:50%;" />

```
library(NMF)

# raw_count2: 行是基因，列是样本

# 过滤掉表达很低的基因（比如在多数样本中都是0）
keep_genes <- rowSums(raw_count2 > 1) >= 2
expr_mat <- raw_count2[keep_genes, ]

#log转换一下
log_expr <- log2(expr_mat + 1)

# 设定重复次数（越高越稳定，但越慢）
nmf_res <- nmf(log_expr, 2:10, nrun = 30, seed = 123)

# 可视化 rank survey 图
plot(nmf_res)

##rank3是最优解 
best_fit <- nmf(log_expr, rank = 2, nrun = 30, seed = 123)

consensusmap(best_fit)
```



感觉还可以走一下deseq2

```
library(DESeq2)
condition_df <- data.frame(condition = factor(condition$condition))
condition_df$condition <- relevel(condition_df$condition, ref = "TREM2_WT-young")
rownames(condition_df) <- rownames(condition)

# 确保 count matrix 列顺序和 condition 行顺序一致
raw_count_filtered <- raw_count2[, rownames(condition_df)]

# 创建 DESeq2 对象
dds <- DESeqDataSetFromMatrix(countData = raw_count_filtered,
                              colData = condition_df,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)


resOrdered <- res[order(res$padj), ]
head(resOrdered)

> head(resOrdered)
log2 fold change (MLE): condition TREM2 KO.young vs TREM2 WT.young 
Wald test p-value: condition TREM2 KO.young vs TREM2 WT.young 
DataFrame with 6 rows and 6 columns
           baseMean log2FoldChange     lfcSE      stat      pvalue        padj
          <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
Trem2      175.3522       -2.58066  0.183979 -14.02695 1.06639e-44 7.15335e-41
Tgfbr1     756.2895        1.00366  0.122466   8.19545 2.49665e-16 8.37376e-13
Sft2d1      89.2022        4.03133  0.496065   8.12662 4.41429e-16 9.87035e-13
Ccr5       330.0379        1.50551  0.188190   7.99994 1.24478e-15 2.08750e-12
Rps12-ps9  176.7929       -1.58750  0.206829  -7.67539 1.64916e-14 2.21252e-11
Ctnnd1      36.6419        2.30624  0.330627   6.97535 3.05115e-12 3.41118e-09
```

做个火山图看看

```
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

<img src="C:\Users\dengy\AppData\Roaming\Typora\typora-user-images\image-20250415101359408.png" alt="image-20250415101359408" style="zoom:50%;" />

还可以 看看取交集会怎么样 

```
DEG_UP <-subset(res,res$log2FoldChange>=1 & res$pvalue <=0.05)
head(DEG_UP)
log2 fold change (MLE): condition TREM2 KO.old vs TREM2 WT.old 
Wald test p-value: condition TREM2 KO.old vs TREM2 WT.old 
DataFrame with 6 rows and 6 columns
        baseMean log2FoldChange     lfcSE      stat      pvalue        padj
       <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
Abracl   24.7946        1.44267  0.326777   4.41484 1.01084e-05 2.16613e-03
Anxa1    11.8257        2.74638  0.680380   4.03653 5.42469e-05 8.18023e-03
Ccr1     32.0591        1.24887  0.313723   3.98079 6.86861e-05 9.81247e-03
Ccr2     56.6624        1.13223  0.252888   4.47719 7.56313e-06 1.71074e-03
Cd209a   23.3102        2.54370  0.421850   6.02988 1.64086e-09 1.48461e-06
Cd209f   44.5710        1.06053  0.272000   3.89902 9.65835e-05 1.15659e-02

#去除非编码基因 
DEG_UP$GeneID <-row.names(DEG_UP)
DEG_UP <-as.data.frame(DEG_UP)
filtered_DEG <- DEG_UP %>%
  dplyr::filter(!grepl("^ens", GeneID, ignore.case = TRUE) & 
                !grepl("rik$", GeneID, ignore.case = TRUE) &
                !grepl("^rps", GeneID, ignore.case = TRUE) &
                !grepl("^rpl", GeneID, ignore.case = TRUE) &
                !grepl("^gm", GeneID, ignore.case = TRUE))
DEG_UP_young <-DEG_UP 
DEG_UP_old <-DEG_UP 
inter_gene <- intersect(row.names(DEG_UP_old), row.names(DEG_UP_young))                



gene = bitr(inter_gene, #数据集
fromType="SYMBOL", #输入为SYMBOL格式
toType="ENTREZID",  # 转为ENTERZID格式
OrgDb="org.Mm.eg.db") #小鼠数据库 

head(gene)
  SYMBOL ENTREZID
2 Abca13   268379
3 Abracl    73112
4 Adam19    11492
5  Adam9    11502
6  Anxa1    16952
7  Apba3    57267


ego_ALL <- enrichGO(gene = gene$ENTREZID,
                OrgDb = org.Mm.eg.db, #没有organism="human"，改为OrgDb=org.Hs.eg.db
				#keytype = 'ENSEMBL',
                ont = "ALL", #也可以是 CC  BP  MF中的一种
                pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                pvalueCutoff = 1, #P值会过滤掉很多，可以全部输出
                qvalueCutoff = 1,
				readable = TRUE) #Gene ID 转成gene Symbol ，易读




DEG_UP_old$GeneID<-row.names(DEG_UP_old)
DEG_UP_old<-as.data.frame(DEG_UP_old)
DEG_UP_old<- DEG_UP_old %>%
  dplyr::filter(!grepl("^ens", GeneID, ignore.case = TRUE) & 
                !grepl("rik$", GeneID, ignore.case = TRUE) &
                !grepl("^rps", GeneID, ignore.case = TRUE) &
                !grepl("^rpl", GeneID, ignore.case = TRUE) &
                !grepl("^rp23", GeneID, ignore.case = TRUE) &
                !grepl("^rp24", GeneID, ignore.case = TRUE) &
                !grepl("^mt-", GeneID, ignore.case = TRUE) &
                !grepl("^gm", GeneID, ignore.case = TRUE))

DEG_UP_young$GeneID<-row.names(DEG_UP_young)
DEG_UP_young<-as.data.frame(DEG_UP_young)
DEG_UP_young<- DEG_UP_young %>%
  dplyr::filter(!grepl("^ens", GeneID, ignore.case = TRUE) & 
                !grepl("rik$", GeneID, ignore.case = TRUE) &
                !grepl("^rps", GeneID, ignore.case = TRUE) &
                !grepl("^rpl", GeneID, ignore.case = TRUE) &
                !grepl("^rp23", GeneID, ignore.case = TRUE) &
                !grepl("^rp24", GeneID, ignore.case = TRUE) &
                !grepl("^mt-", GeneID, ignore.case = TRUE) &
                !grepl("^gm", GeneID, ignore.case = TRUE))


inter_gene <- intersect(row.names(DEG_UP_old), row.names(DEG_UP_young))      


venn_list <- list(Trem2KO_old = rownames(DEG_UP_old),
                  Trem2KO_young = rownames(DEG_UP_young))

ggVennDiagram(venn_list, category.names = c("1", "2"), # 设定样本名称
              label = "both", # 可选："both", "count", "percent", "none"
              label_color = "black",
              label_alpha = 0, # 去除文字标签底色
              #edge_lty = "dashed", # 圆圈线条虚线
              edge_size = 1)+ scale_fill_gradient(low = "white", high = alpha("#f8766d", 0.9), name = "gene count")


gene = bitr(row.names(DEG_UP_old), #数据集
fromType="SYMBOL", #输入为SYMBOL格式
toType="ENTREZID",  # 转为ENTERZID格式
OrgDb="org.Mm.eg.db") #小鼠数据库 
ego_ALL <- enrichGO(gene = gene$ENTREZID,
                OrgDb = org.Mm.eg.db, #没有organism="human"，改为OrgDb=org.Hs.eg.db
				#keytype = 'ENSEMBL',
                ont = "ALL", #也可以是 CC  BP  MF中的一种
                pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                pvalueCutoff = 1, #P值会过滤掉很多，可以全部输出
                qvalueCutoff = 1,
				readable = TRUE) #Gene ID 转成gene Symbol ，易读

# 加载必要的R包
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("org.Mm.eg.db", "AnnotationDbi"))

library(org.Mm.eg.db)
library(AnnotationDbi)

go_id <- c("GO:0006935","GO:0010573",
"GO:0001935"
)

symbols <- get(go_id, org.Mm.egGO2ALLEGS) %>% mget(org.Mm.egSYMBOL) %>% unlist()




DEG_UP_young<-as.data.frame(res)
DEG_UP_young$GeneID<-row.names(DEG_UP_young)
DEG_UP_young<- DEG_UP_young %>%
  dplyr::filter(!grepl("^ens", GeneID, ignore.case = TRUE) & 
                !grepl("rik$", GeneID, ignore.case = TRUE) &
                !grepl("^rps", GeneID, ignore.case = TRUE) &
                !grepl("^rpl", GeneID, ignore.case = TRUE) &
                !grepl("^rp23", GeneID, ignore.case = TRUE) &
                !grepl("^rp24", GeneID, ignore.case = TRUE) &
                !grepl("^mt-", GeneID, ignore.case = TRUE) &
                !grepl("^gm", GeneID, ignore.case = TRUE))


DEG_UP_young$status <- ifelse(DEG_UP_young$log2FoldChange >= 1 & DEG_UP_young$pvalue <= 0.05, "Up",
                          ifelse(DEG_UP_young$log2FoldChange <= -1 & DEG_UP_young$pvalue <= 0.05, "Down", "NotSig"))



symbols<-as.data.frame(symbols)
DEG_UP_young %>% inner_join(symbols %>% select(SYMBOL, ENTREZID), by = c("GeneID" = "SYMBOL"))

DEG_UP_young <- DEG_UP_young %>%
  mutate(status = case_when(
    log2FoldChange >= 1 & pvalue <= 0.05 ~ "Up",
    log2FoldChange <= -1 & pvalue <= 0.05 ~ "Down",
    GeneID %in% symbols$SYMBOL ~ "GO: cell chemotaxis",
    TRUE ~ "NotSig"
  ))




for (i in 1:nrow(DEG_UP_young)) {
  gene_id <- DEG_UP_young$GeneID[i]
  
  # 检查 GeneID 是否在 symbols$SYMBOL 中
  if (gene_id %in% symbols) {
    DEG_UP_young$status[i] <- "GO: cell chemotaxis"
  } else {
    # 根据 log2FoldChange 和 pvalue 更新 status
    if (!is.na(DEG_UP_young$log2FoldChange[i]) && !is.na(DEG_UP_young$pvalue[i])) {
      if (DEG_UP_young$log2FoldChange[i] >= 1 & DEG_UP_young$pvalue[i] <= 0.05) {
        DEG_UP_young$status[i] <- "Up"
      } else if (DEG_UP_young$log2FoldChange[i] <= -1 & DEG_UP_young$pvalue[i] <= 0.05) {
        DEG_UP_young$status[i] <- "Down"
      } else {
        DEG_UP_young$status[i] <- "NotSig"
      }
    } else {
      DEG_UP_young$status[i] <- "NotSig"
    }
  }
}





























label_data <- subset(df, status == "GO: cell chemotaxis" & abs(log2FoldChange) >= 1 & pvalue <= 0.05 )
```

