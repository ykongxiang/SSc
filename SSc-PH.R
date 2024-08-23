library(Seurat)
library(data.table)
library(tibble)
library(dplyr)
library(scater)
library(scran)
library(cowplot)
library(rhdf5)
library(Biobase)
library(ggplot2)
load("H:/bioinformatics/transcriptomics/SSc scRNA-RNAseq/SSc数据/GSE244781_scRNA_SSc_ILD/seurat_list.RData")
# 设置工作目录
setwd('H:/bioinformatics/transcriptomics/SSc scRNA-RNAseq/SSc数据/GSE244781_scRNA_SSc_ILD')

# 列出所有 h5 文件
h5_files <- list.files(pattern = "\\.h5$", full.names = TRUE)

# 从文件名中提取样本名称（假设文件名格式为 "GSMXXXXXX_*.h5"）
sample_names <- sub("\\./", "", sub("_.*", "", h5_files))

# 创建一个包含所有文件的名称与其对应的 GSM 的数据框
file_gsm_mapping <- data.frame(
  SampleName = sample_names,
  GSM = sub("\\./", "", sub("_.*", "", basename(h5_files)))
)

# 根据 group_data 中的 GSM 筛选出需要的样本
selected_samples <- group_data$GSM
selected_files <- file_gsm_mapping$SampleName[file_gsm_mapping$GSM %in% selected_samples]

# 根据筛选出的样本名筛选 h5 文件
selected_h5_files <- h5_files[file_gsm_mapping$SampleName %in% selected_files]

# 初始化 Seurat 对象列表
seurat_list <- list()

# 循环读取每个 h5 文件并创建 Seurat 对象
for (i in 1:length(selected_h5_files)) {
  # 读取数据
  scRNA_data <- Read10X_h5(file = selected_h5_files[i])
  
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
                                   project = basename(selected_h5_files[i]),
                                   min.features = 200,
                                   min.cells = 3)
  
  # 将处理后的 Seurat 对象添加到列表
  seurat_list[[i]] <- seurat_obj
  
  # 清理内存
  gc()
}

# 继续进行数据整合和分析（如前所述）
# 初始化合并对象
SSc_PH <- seurat_list[[1]]
# 循环合并剩余的 Seurat 对象
for (i in 2:length(seurat_list)) {
  SSc_PH <- merge(x = SSc_PH, 
                        y = seurat_list[[i]], 
                        add.cell.ids = names(seurat_list)[i], 
                        project = "SSc_combined_project")
}
rm(scRNA_data,counts_matrix,seurat_list,seurat_obj)
save(SSc_PH,file='SSc_PH.RData')
#anchor方法 ------------------------------------------------------------------------------
# 找到整合的锚点
save(SSc_PH,file = 'SSc_PH.RData')
anchors <- FindIntegrationAnchors(object.list = SSc_PH, dims = 1:30)
save(anchors,file = 'anchors.RData')
# 整合数据
SSc_PH <- IntegrateData(anchorset = anchors, dims = 1:30)
DimPlot(combined,dims= 1:30)
#harmony----------------------------------------------------------------------------
library(harmony)
library(Seurat)
load("H:/bioinformatics/transcriptomics/SSc scRNA-RNAseq/SSc数据/GSE244781_scRNA_SSc_ILD/seurat_list.RData")
save(SSc_PH,file = 'SSc_PH.RData')
SSc_PH <- RunHarmony(
  object = SSc_PH,
  group.by.vars = "",  # 或者是你的分组变量名
  dims.use = 1:20,               # 使用的PCA维度
  assay.use = "RNA",             # 使用的assay层
  plot_convergence = TRUE        # 画出收敛图像
)
#------------------------------------------------------------------------------
load('H:/bioinformatics/transcriptomics/SSc scRNA-RNAseq/SSc数据/GSE244781_scRNA_SSc_ILD/SSc_PH.RData')
group_data <- data.frame(
  GSM = c("GSM7829642", "GSM7829641", "GSM7829640", "GSM7829639", "GSM7829638",
          "GSM7829637", "GSM7829636", "GSM7829635", "GSM7829634", "GSM7829633",
          "GSM7829632",  "GSM7829650", "GSM7829649", "GSM7829648", "GSM7829647",
          "GSM7829646", "GSM7829645", "GSM7829644", "GSM7829643"),
  Sex = c("Male", "Male", "Male", "Male", "Male", "Male", "Female", "Female", "Female", "Female",
          "Male",  "Female","Female", "Female", "Female", "Female", "Male", "Female", "Female"),
  Disease = c("SSc", "SSc", "SSc", "SSc", "SSc", "SSc", "SSc", "SSc", "iPAH", "iPAH",'iPAH', "SSc", "SSc",
              "SSc", "SSc", "SSc", "SSc", "SSc", "SSc")
)
SSc_PH@meta.data[["orig.ident"]] <- sub("_.*", "", SSc_PH@meta.data[["orig.ident"]])
SSc_PH$Sex <- group_data$Sex[match(SSc_PH$orig.ident, group_data$GSM)]
SSc_PH$Disease <- group_data$Disease[match(SSc_PH$orig.ident, group_data$GSM)]

SSc_PH[["percent.mt"]] <- PercentageFeatureSet(SSc_PH, pattern = "^MT-")
fivenum(SSc_PH@meta.data[["percent.mt"]])
VlnPlot(SSc_PH,features = "percent.mt",pt.size = 0,group.by = "orig.ident") 
SSc_PH <- subset(SSc_PH, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
VlnPlot(SSc_PH,features = "percent.mt",pt.size = 0,group.by = "orig.ident")
SSc_PH <- NormalizeData(SSc_PH)
SSc_PH <- FindVariableFeatures(SSc_PH, selection.method = "vst", nfeatures = 2000)
SSc_PH <- ScaleData(SSc_PH)
SSc_PH <- RunPCA(SSc_PH)
SSc_PH <- FindNeighbors(SSc_PH, dims = 1:30, reduction = "pca")
SSc_PH <- FindClusters(SSc_PH, resolution = 2, cluster.name = "unintegrated_clusters")
SSc_PH <- RunUMAP(SSc_PH, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# 筛选女性样本
SSc_PH_Female <- subset(SSc_PH, subset = Sex == "Female")
SSc_PH_Female <- RunPCA(SSc_PH_Female)
SSc_PH_Female <- RunUMAP(SSc_PH_Female, dims = 1:30)
DimPlot(SSc_PH_Female, reduction = "umap", group.by = "Disease", label = TRUE, repel = TRUE)

DimPlot(SSc_PH,reduction = "umap.unintegrated",group.by = 'Disease',split.by = 'Sex',label = TRUE)
DimPlot(SSc_PH, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters"))
SSc_PH <- IntegrateLayers(object = SSc_PH,method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                           verbose = FALSE)
SSc_PH <- SCTransform(SSc_PH, verbose = FALSE)
SSc_PH <- IntegrateLayers(SSc_PH, method = SCTIntegration, dims = 1:30, orig.reduction = "pca", new.reduction = "integrated.sct")
SSc_PH[["RNA"]] <- JoinLayers(SSc_PH[["RNA"]])
SSc_PH <- FindNeighbors(SSc_PH, dims = 1:30, reduction = "integrated.cca")
SSc_PH <- FindClusters(SSc_PH, resolution = 1)
SSc_PH <- RunUMAP(SSc_PH, dims = 1:30, reduction = "integrated.cca")
DimPlot(SSc_PH, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"))
DimPlot(SSc_PH, reduction = "umap", split.by = "orig.ident")
save(SSc_PH,file = 'SSc_PH.RData')
#-------------------------------------------------------------------------------------------
SSc_PH.markers <- FindAllMarkers(SSc_PH,only.pos = T ,min.pct = 0.1 ,logfc.threshold = 0.25)
top10 <- SSc_PH.markers%>%group_by(cluster)%>%top_n(n=10, wt=avg_log2FC)
DoHeatmap(SSc_PH,features = top10$gene)
marker_genes <- list(
  Tip = c('COL4A1', 'HSPG2', 'VWA1'),
  Phalanx  = c('CXCL12', 'COL15A1'),
  Arterial = c('IGFBP3', 'DKK2', 'GJA5'),
  Systemic_Venous = c('CCL14', 'ACKR1', 'SELE'),
  Pulmonary_Venous = c('CCL23', 'ACKR1', 'CPE'),
  Aerocyte = c('HPGD', 'IL1RL1'),#'SOSTDC1'
  General_Capillary = c('FCN3', 'IL7R', 'TMEM100'),
  Proliferating_Vascular = c('PCLAF', 'STMN1', 'PTTG1'),
  Transitional = c('NDUFA4L2', 'AGT','COX4I1'),#'COX4I1'
  Lymphatic = c('CCL21', 'TFF3', 'MMRN1'),
  Pericyte = c('COX4I1', 'NDUFA4L2', 'HIGD1B'),
  Smooth_Muscle = c('ACTA2', 'TAGLN', 'TPM2'),
  Adventital_Fibroblast = c('PLA2G2A', 'SFRP2', 'CFD'),
  Alveolar_Fibroblast = c('G0S2', 'RARRES2', 'GPC3'),
  WIF1_Fibroblast = c('APOD', 'APOE', 'ASPN'),
  Myofibroblast = c('COL1A1', 'COL3A1', 'COL1A2')
)
all_marker_genes <- unlist(marker_genes)
DefaultAssay(SSc_PH) <- "RNA"

# 使用 DotPlot 进行可视化
DotPlot(SSc_PH, features = unique_marker_genes, group.by = 'seurat_clusters', dot.scale = 6) +
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = round(avg.exp, 2)), size = 2)

unique_marker_genes <- unique(all_marker_genes)

# 绘制热图
DoHeatmap(SSc_PH, features = unique_marker_genes, group.by = 'seurat_clusters') 
cell_type <- c('0,1'='WIF1_Fibroblast',
               '2,3,4'='',
               '11,25'= 'Systemic_Venous',
               '10,20'= 'Myofibroblast',
               '')
cell_names <- colnames(SSc_PH)

# 提取样本号（orig.ident信息）
sample_ids <- SSc_PH$orig.ident

# 创建新的细胞名称，将样本号添加到细胞名称前
new_cell_names <- paste(sample_ids, cell_names, sep = "_")

# 更新 Seurat 对象中的细胞名称
colnames(SSc_PH) <- new_cell_names

library(ggplot2)
library(reshape2)

library(Seurat)
library(ggplot2)
library(reshape2)

plot_marker_heatmap <- function(seurat_obj, marker_genes, colors = c("blue", "white", "red"), facet_nrow = 1) {
  
  # 定义检测函数
  check_marker_genes <- function(seurat_obj, marker_genes) {
    # 提取 Seurat 对象中的所有基因名称
    all_genes <- rownames(seurat_obj) # 使用 Seurat 对象中的 RNA counts 矩阵
    
    # 扁平化 marker_genes 列表
    all_marker_genes <- unlist(marker_genes)
    
    # 检查标记基因是否在 Seurat 对象中
    missing_genes <- setdiff(all_marker_genes, all_genes)
    
    # 打印结果并处理缺失基因
    if (length(missing_genes) > 0) {
      cat("Missing genes:\n")
      print(missing_genes)
      stop("Error: Some marker genes are missing from the Seurat object. Please check the gene names.")
    } else {
      cat("All marker genes are present in the Seurat object.\n")
    }
  }
  
  # 检查标记基因是否存在
  check_marker_genes(seurat_obj, marker_genes)
  
  marker_genes_df <- data.frame()
  
  # 遍历列表中的每个元素，并将其添加到数据框中
  for (cell_type in names(marker_genes)) {
    genes <- marker_genes[[cell_type]]
    temp_df <- data.frame(CellType = cell_type, GeneName = genes)
    marker_genes_df <- rbind(marker_genes_df, temp_df)
  }
  
  # 提取标记基因表达数据并转换为 data.frame
  expr_data <- FetchData(SSc_PH, vars = marker_genes_df$GeneName)
  expr_data <- as.data.frame(expr_data)
  expr_data$cell <- rownames(expr_data)
  
  # 提取细胞的簇信息并确保是标准的因子格式
  clusters <- Idents(SSc_PH)
  expr_data$cluster <- as.factor(clusters[rownames(expr_data)])
  # 将数据转换为长格式

  # 创建热图
  p <- ggplot(expr_data_long, aes(x = cell, y = gene, fill = expression)) +
    geom_tile() +
    scale_fill_gradientn(colors = colors) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),  # 隐藏 x 轴标签以增强可读性
          axis.ticks.x = element_blank(),
          panel.spacing = unit(0.1, "lines")) +
    labs(x = "Cells", y = "Genes", fill = "Expression") +
    facet_wrap(~cluster, scales = "free_x", nrow = facet_nrow) +  # 使用 facet_wrap 确保每个簇等宽展示
    scale_x_discrete(expand = expansion(mult = c(0.05, 0.05)))  # 调整每个簇的宽度
  
  return(p)
}
write.csv(expr_data_long,file = 'expr_data_long.csv')

p <- plot_marker_heatmap(SSc_PH, marker_genes = marker_genes, colors = c("blue", "white", "red"), facet_nrow = 1)

