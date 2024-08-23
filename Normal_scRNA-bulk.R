library(data.table)
library(dplyr)     
library(ggplot2)   
library(pheatmap)   
library(DESeq2)    
library(edgeR)     
library(limma)     
library(tinyarray) 
library(org.Hs.eg.db)
library(readr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(KEGGREST)
library(rlang)
library(KEGGgraph)
library(pheatmap)
run_limma_analysis <- function(data, design, contrast_vector, topNum=Inf) {
  
  # 进行线性建模
  fit <- lmFit(data, design)
  
  # 进行统计检验
  fit <- eBayes(fit)
  
  # 创建对比矩阵并计算对比
  contrast.matrix <- makeContrasts(contrasts = contrast_vector, levels = design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit)
  
  # 获取差异表达基因结果
  result_limma <- na.omit(topTable(fit, number = topNum))
  
  return(result_limma)
}
plot_volcano <- function(result, log2FoldChange = "logFC", pvalue = "adj.P.Val", rank = 10) {
  # 添加显著性分类
  result$significance <- "Not significant"
  result$significance[result[[log2FoldChange]] > 1 & result[[pvalue]] < 0.05] <- "Upregulated"
  result$significance[result[[log2FoldChange]] < -1 & result[[pvalue]] < 0.05] <- "Downregulated"
  
  # 筛选出显著的基因
  significant_genes <- result[result$significance != "Not significant", ]
  
  # 按差异度排序并选择前 rank 个基因
  top_genes <- significant_genes[order(-abs(significant_genes[[log2FoldChange]])), ]
  top_genes <- head(top_genes, rank)
  
  # 绘制火山图
  volcano_plot <- ggplot(result, aes(x = !!sym(log2FoldChange), y = -log10(!!sym(pvalue)))) +
    geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Downregulated" = "blue2", 
                                  "Not significant" = "gray30", 
                                  "Upregulated" = "red2")) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          legend.position = c(0.85, 0.1)) +
    labs(x = 'log2 Fold Change', y = '-log10 Adjusted P-value', color = 'Significance') +
    geom_vline(xintercept = c(-1, 1), color = 'gray', linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05), color = 'gray', linewidth = 0.5) +
    geom_text_repel(data = top_genes, aes(label = rownames(top_genes)), 
                    size = 3, box.padding = 0.3, point.padding = 0.2)
  
  # 显示火山图
  print(volcano_plot)
}

# 定义样本ID和对应的分组信息
sample_ids <- c("GSM1003068", "GSM1003069", "GSM1003070", "GSM1003071",
                "GSM1003072", "GSM1003073", "GSM1003074", "GSM1003075",
                "GSM1003076", "GSM1003077", "GSM1003078",
                "GSM1003058", "GSM1003067", "GSM1003059", "GSM1003060",
                "GSM1003061", "GSM1003062", "GSM1003063", "GSM1003064",
                "GSM1003065", "GSM1003066")

group_assignments <- c(rep("SSc_ILD", 8), rep("UIP", 3), rep("Control", 10))

# 创建分组信息数据框
group_info <- data.frame(
  Sample_ID = sample_ids,
  Group = factor(group_assignments, levels = c("SSc_ILD", "UIP", "Control"))
)

GDS <- read.table("H:/bioinformatics/transcriptomics/SSc scRNA-RNAseq\\GDS4995_full.soft",sep = "\t", comment.char = "#", stringsAsFactors = FALSE, 
                  header = FALSE, fill = TRUE)
Exp <- GDS[,c(3:23,25)]
Exp <- Exp[-1,]
Exp1 <- sapply(Exp[,-22], function(x) as.numeric(as.character(x)))
colnames(Exp) <- c(sample_ids,'Symbol')
Exp <- aggregate(Exp1, by = list(Exp$Symbol), FUN = mean) 
rownames(Exp) <- Exp[,1]
Exp1 <- na.omit(Exp)

group <- group_info$Group
DGE_exp <- DGEList(counts = Exp1, group = group)
design <- model.matrix(~0 + group, data = DGE_exp$samples)
colnames(design) <- levels(group)

exp_brca_norm <- log2(DGE_exp$counts + 1)
keep.exprs <- filterByExpr(DGE_exp)
DGE_exp <- DGE_exp[keep.exprs, , keep.lib.sizes = FALSE]
exp_brca_norm_filtered <- log2(DGE_exp$counts + 1)

# 绘制过滤后数据分布图
plot(density(exp_brca_norm_filtered[,1]), col = col[1], lwd = 2, ylim = c(0, 0.26), las = 2, main = "Filtered Data", xlab = "Log2 CPM")
abline(v = lcpm.cutoff, lty = 3, col = "red")
for (i in 2:nsamples) {
  den <- density(exp_brca_norm_filtered[,i])
  lines(den$x, den$y, col = col[i], lwd = 2)
}
legend("topright", legend = c("Sample 1", paste("Sample", 2:nsamples)), text.col = col, bty = "n", cex = 0.8)
par(mfrow = c(1, 1))

par(mfrow = c(2, 3))
boxplot(DGE_exp$counts, outline = FALSE, col = col, main = "A. Unnormalised", ylab = "Raw Count")
norm_factors <- calcNormFactors(DGE_exp, method = "TMM")
boxplot(cpm(DGE_exp$counts, normalized.lib.sizes = norm_factors$lib.size), outline = FALSE, col = col, main = "B. TMM Normalised", ylab = "Normalised Count")
cpm_counts <- cpm(DGE_exp$counts, normalized.lib.sizes = norm_factors$lib.size)
boxplot(cpm_counts, outline = FALSE, col = col, main = "C. CPM", ylab = "CPM")
log_cpm_counts <- cpm(DGE_exp$counts, log = TRUE, normalized.lib.sizes = norm_factors$lib.size)
boxplot(log_cpm_counts, outline = FALSE, col = col, main = "D. Log-CPM", ylab = "Log-CPM")
v <- voom(DGE_exp, design, plot = FALSE)
voom_counts <- v$E
boxplot(voom_counts, outline = FALSE, col = col, main = "E. Voom Normalised Log-CPM", ylab = "Log-CPM")

par(mfrow = c(1, 1))
Exp3 <- voom_counts
colnames(Exp3) <- sample_ids
plotMDS(Exp3, col = rep(c('red', 'blue', 'black'), each = 3))

DGE_exp <- DGEList(counts = Exp3, group = group)

exp_brca_norm <- log2(Exp3 + 1) 
rownames(exp_brca_norm) <- rownames(Exp3)

result_UIP <- run_limma_analysis(exp_brca_norm, design, "Control - UIP")
result_SSc_ILD <- run_limma_analysis(exp_brca_norm, design, "Control - SSc_ILD")
result_UIP_SSc <- run_limma_analysis(exp_brca_norm, design, "UIP - SSc_ILD")
plot_volcano(result_UIP, rank = 10)
plot_volcano(result_SSc_ILD, rank = 10)
plot_volcano(result_UIP_SSc, rank = 10)

#------------------------------------------
# 整合的富集分析函数
perform_enrichment_analysis <- function(diff_genes, 
                                        fromType = "SYMBOL", 
                                        toType = "ENTREZID",
                                        orgDb = org.Hs.eg.db, 
                                        pvalueCutoff = 0.05, 
                                        qvalueCutoff = 0.2, 
                                        ont = "BP") {
  
  # 1. 基因符号转换为ENTREZ ID
  convert_gene_symbols_to_entrez <- function(diff_genes_symbols) {
    # 清理包含 "///" 的基因符号
    diff_genes_symbols_clean <- diff_genes_symbols %>%
      mutate(Gene.Symbol = sapply(strsplit(Gene.Symbol, "///"), `[`, 1))
    
    # 转换基因符号为 ENTREZ ID
    gene_conversion <- bitr(diff_genes_symbols_clean$Gene.Symbol, 
                            fromType = fromType, 
                            toType = toType, 
                            OrgDb = orgDb)
    
    # 合并转换结果
    diff_genes_entrez <- merge(diff_genes_symbols_clean, gene_conversion, 
                               by.x = "Gene.Symbol", by.y = fromType)
    
    return(diff_genes_entrez$ENTREZID)
  }
  
  # 提取差异基因的符号名称
  diff_genes_symbols <- data.frame(Gene.Symbol = rownames(diff_genes), stringsAsFactors = FALSE)
  
  # 转换基因符号为ENTREZ ID
  entrez_ids <- convert_gene_symbols_to_entrez(diff_genes_symbols)
  
  # 2. GO 富集分析
  perform_go_enrichment <- function(entrez_ids) {
    go_enrich <- enrichGO(gene = entrez_ids, 
                          OrgDb = orgDb, 
                          keyType = toType, 
                          ont = ont,  # BP: Biological Process
                          pvalueCutoff = pvalueCutoff, 
                          qvalueCutoff = qvalueCutoff)
    
    return(go_enrich)
  }
  
  # 3. KEGG 富集分析
  perform_kegg_enrichment <- function(entrez_ids) {
    kegg_enrich <- enrichKEGG(gene = entrez_ids,
                              organism = "hsa",
                              pvalueCutoff = pvalueCutoff,
                              qvalueCutoff = qvalueCutoff)
    
    return(kegg_enrich)
  }
  
  # 运行GO和KEGG富集分析
  go_results <- perform_go_enrichment(entrez_ids)
  kegg_results <- perform_kegg_enrichment(entrez_ids)
  
  # 返回结果
  return(list(GO = go_results, KEGG = kegg_results))
}


enrichment_SSc_ILD <- perform_enrichment_analysis(result_SSc_ILD)
enrichment_UIP_SSc <- perform_enrichment_analysis(result_UIP_SSc)
enrichment_UIP <- perform_enrichment_analysis(result_UIP)
dotplot(enrichment_UIP_SSc$GO, showCategory = 20)
dotplot(enrichment_UIP_SSc$KEGG, showCategory = 20)
# 查看GO和KEGG富集分析的结果
head(enrichment_results$GO)
head(enrichment_results$KEGG)
GO_SSc_ILD <- as.data.frame(enrichment_results$GO)


#scRNA-seq-----------------------------------------
library(Seurat)
if (!require(tidyverse)) install.packages("tidyverse")
library(ggplot2)
library(SingleR)
library(GEOquery)
library(Matrix)
library(hdf5r)
library(tools)
library(harmony)
library(glmGamPoi)
library(celldex)
data <- data.frame(
  AGE = c(76, 56, 56, 55, 57, 18, 23, 23, 21, 33, 64, 38, 62, 63, 81, 21, 50, 36, 67, 67, 43, 43, 53, 53, 64, 64,64, 43,33, 62, 59,42, 58, 54),
  `Library Name` = c("GSM7829618", "GSM7829619", "GSM7829619", "GSM7829620", "GSM7829621", "GSM7829622", 
                     "GSM7829623", "GSM7829624", "GSM7829625", "GSM7829626", "GSM7829627", "GSM7829628",
                     "GSM7829629", "GSM7829630", "GSM7829631", "GSM7829632", "GSM7829633", "GSM7829634",
                     "GSM7829635", "GSM7829636", "GSM7829637", "GSM7829638", "GSM7829639", "GSM7849640",
                     "GSM7849641", "GSM7849642", "GSM7849643", "GSM7849644", "GSM7849645", "GSM7849646",
                     "GSM7849647", "GSM7849648", "GSM7849649", "GSM7849650"),
  `disease_state` = c("normal control", "normal control", "normal control", "normal control", "normal control",
                      "normal control", "normal control", "normal control", "normal control", "normal control",
                      "normal control", "normal control", "systemic sclerosis", "systemic sclerosis",
                      "systemic sclerosis", "systemic sclerosis", "systemic sclerosis", "systemic sclerosis",
                      "systemic sclerosis", "systemic sclerosis", "idiopathic pulmonary arterial hypertension",
                      "idiopathic pulmonary arterial hypertension", "idiopathic pulmonary arterial hypertension",
                      "systemic sclerosis", "systemic sclerosis", "systemic sclerosis", "systemic sclerosis",
                      "systemic sclerosis", "systemic sclerosis", "systemic sclerosis", "systemic sclerosis",
                      "systemic sclerosis", "systemic sclerosis"),
  Sex = c("male", "male", "male", "male", "female", "male", "female", "female", "male", "male", "female",
          "female", "male", "male", "male", "female", "female", "female", "female", "male", "male",
          "male", "male", "male", "male", "female", "female", "male", "female", "female", "female",
          "female", "female"),
  `lung_section` = c(0, 0, 0, 0, 0, 0,"Lower lobe", "Upper lobe", 0, 0, 0, "50% CD45, EPCAM depleted","50% CD45, EPCAM depleted",
                     0, 0, 0, 0,0, "Lower lobe",
                     "Upper lobe", "Lower lobe", "Upper lobe", "Lower lobe", "Upper lobe", "Lower lobe",
                     "Upper lobe", 0, 0, 0,0,0, "50% CD45, EPCAM depleted", "50% CD45 depleted", "50% CD163 depleted")
)

setwd('H:/bioinformatics/transcriptomics/SSc scRNA-RNAseq/SSc数据/GSE244781_scRNA_SSc_ILD')
h5_files <- list.files(pattern = "\\.h5$", full.names = TRUE)
sample_names <- sub("\\./", "", sub("_.*", "", h5_files))
seurat_list <- list()
library(Seurat)

# 设置工作目录
setwd('H:/bioinformatics/transcriptomics/SSc scRNA-RNAseq/SSc数据/GSE244781_scRNA_SSc_ILD')

# 获取所有 .h5 文件
h5_files <- list.files(pattern = "\\.h5$", full.names = TRUE)
sample_names <- sub("\\./", "", sub("_.*", "", h5_files))

# 初始化 Seurat 对象列表
seurat_list <- list()

# 循环读取每个 h5 文件
for (i in 1:length(h5_files)) {
  # 读取数据
  scRNA_data <- Read10X_h5(file = h5_files[i])
  
  # 确定数据层
  if (is.list(scRNA_data)) {
    if ("Gene Expression" %in% names(scRNA_data)) {
      counts_matrix <- scRNA_data[["Gene Expression"]]
    } else {
      counts_matrix <- scRNA_data[[1]]
    }
  } else {
    counts_matrix <- scRNA_data
  }
  
  # 创建 Seurat 对象
  seurat_obj <- CreateSeuratObject(counts = counts_matrix, 
                                   project = sample_names[i],
                                   min.features = 200,
                                   min.cells = 3)
  
  # 数据预处理
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  
  # 将处理后的 Seurat 对象添加到列表
  seurat_list[[i]] <- seurat_obj
  
  # 清理内存
  gc()
}

#anchor方法 ------------------------------------------------------------------------------
# 找到整合的锚点
anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:50)
combined <- IntegrateData(anchorset = anchors, dims = 1:30)
DimPlot(combined,dims= 1:50)
#harmony----------------------------------------------------------------------------
library(harmony)
library(Seurat)

save(combined,file = 'combined.RData')
combined_filtered <- subset(combined, cells = WhichCells(combined, idents = group_data[[1]]))
rm(combined)
combined <- RunHarmony(
  object = seurat_list,
  group.by.vars = "orig.ident",  # 或者是你的分组变量名
  dims.use = 1:20,               # 使用的PCA维度
  assay.use = "RNA",             # 使用的assay层
  plot_convergence = TRUE        # 画出收敛图像
)
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:20)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)
save(combined,file = 'seurat_list.RData')
DimPlot(combined, reduction = "umap", group.by = "orig.ident")

combined_filtered <- subset(combined, cells = selected_cells)
save(combined_filtered,file = 'combined_filtered.RData')
combined_filtered <- ScaleData(combined_filtered)
combined_filtered <- RunPCA(combined_filtered, npcs = 20)
combined_filtered <- RunUMAP(combined_filtered, reduction = "harmony", dims = 1:20)
combined_filtered <- FindNeighbors(combined_filtered, reduction = "harmony", dims = 1:20)
combined_filtered <- FindClusters(combined_filtered, resolution = 0.5)
DimPlot(combined_filtered, reduction = "umap", group.by = "orig.ident")
#------------------------------------------------------------------------------
cell_names <- colnames(combined_filtered)

# 将 GSM 提取出来作为样本名
cell_samples <- sub("_.*", "", cell_names)

rownames(group_data) <- group_data$GSM

# 匹配细胞样本名和 `group_data`
cell_annotations <- group_data[cell_samples, ]

# 将 annotations 添加到 Seurat 对象的 metadata 中
combined_filtered$group <- cell_annotations$Sex
combined_filtered$disease <- cell_annotations$Disease
process_seurat <- function(seurat_obj, npcs = 20, dims = 1:20, resolution = 0.5) {
  # 数据标准化
  seurat_obj <- ScaleData(seurat_obj)
  
  # PCA降维
  seurat_obj <- RunPCA(seurat_obj, npcs = npcs)
  
  # UMAP降维
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = dims)
  
  # 查找邻居
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = dims)
  
  # 聚类
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  
  return(seurat_obj)
}
# 按性别子集分析
SSc_male <- subset(combined_filtered, subset = group == "Male" & disease == "SSc")
SSc_female <- subset(combined_filtered, subset = group == "Female" & disease == "SSc")
iPAH_male <- subset(combined_filtered, subset = group == "Male" & disease == "iPAH")
iPAH_female <- subset(combined_filtered, subset = group == "Female" & disease == "iPAH")
SSc_female <- combined_SSc_female
SSc_male <- combined_SSc_male
iPAH_female <- combined_iPAH_female
iPAH_male <- combined_iPAH_male
save(iPAH_female,iPAH_male,SSc_female,SSc_male,file = 'grouped_SSc_PH')
load('grouped_SSc_PH')
# 对每个子集进行预处理
SSc_male <- process_seurat(SSc_male)
SSc_female <- process_seurat(SSc_female)
iPAH_male <- process_seurat(iPAH_male)
iPAH_female <- process_seurat(iPAH_female)

# 进行差异基因分析
markers <- FindMarkers(SSc_female, ident.1 = "iPAH_female", ident.2 = "SSc_female")
