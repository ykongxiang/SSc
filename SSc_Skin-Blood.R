if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(Seurat)
library(data.table)
library(tibble)
library(dplyr)
library(Azimuth)
library(harmony)
setwd('H:/bioinformatics/transcriptomics/SSc scRNA-RNAseq/SSc数据/GSE195452_scRNA_skin_blood')
# 读取所有数据
read_all_data <- function(data_dir) {
  all_data <- Read10X(data.dir = data_dir)
  
  # 提取基因表达数据层
  counts <- all_data[["Gene Expression"]]
  
  return(counts)
}
new_files <- c('GSM5884742-Nor','GSM5884740-SSc')
# 读取第一个样本的数据
counts <- read_all_data(new_files[1])
GSM742 <- CreateSeuratObject(counts = counts,
                             project = 'GSM5884742-Nor',
                             min.features = 200, 
                             min.cells = 3)

# 读取第二个样本的数据
counts <- read_all_data(new_files[2])
GSM740 <- CreateSeuratObject(counts = counts,
                             project = 'GSM5884740-SSc',
                             min.features = 200, 
                             min.cells = 3)
rm(counts)
SSc_Nor <- merge(x = GSM742, y = GSM740, add.cell.ids = c("Nor", "SSc"), project = "SSc_Nor_Merged")
SSc_Nor[["percent.mt"]] <- PercentageFeatureSet(SSc_Nor, pattern = "^MT-")
fivenum(SSc_Nor@meta.data[["percent.mt"]])
VlnPlot(SSc_Nor,features = "percent.mt",pt.size = 0,group.by = "orig.ident") 
SSc_Nor <- subset(SSc_Nor, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
VlnPlot(SSc_Nor,features = "percent.mt",pt.size = 0,group.by = "orig.ident")
SSc_Nor <- NormalizeData(SSc_Nor)
SSc_Nor <- FindVariableFeatures(SSc_Nor, selection.method = "vst", nfeatures = 2000)
SSc_Nor <- ScaleData(SSc_Nor)
SSc_Nor <- RunPCA(SSc_Nor)
SSc_Nor <- FindNeighbors(SSc_Nor, dims = 1:30, reduction = "pca")
SSc_Nor <- FindClusters(SSc_Nor, resolution = 2, cluster.name = "unintegrated_clusters")
SSc_Nor <- RunUMAP(SSc_Nor, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(SSc_Nor, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters"))
SSc_Nor <- IntegrateLayers(object = SSc_Nor,method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE)
SSc_Nor[["RNA"]] <- JoinLayers(SSc_Nor[["RNA"]])
SSc_Nor <- FindNeighbors(SSc_Nor, dims = 1:30, reduction = "integrated.cca")
SSc_Nor <- FindClusters(SSc_Nor, resolution = 1)
SSc_Nor <- RunUMAP(SSc_Nor, dims = 1:30, reduction = "integrated.cca")
DimPlot(SSc_Nor, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"))
DimPlot(SSc_Nor, reduction = "umap", group.by = "seurat_clusters",label = TRUE,label.size = 5)
DimPlot(SSc_Nor, reduction = "umap", split.by = "orig.ident")
save(SSc_Nor,file = 'SSc_Nor.RData')
load('SSc_Nor.RData')

#-------------------------------------------------------------------------------
SSc_Nor.markers <- FindAllMarkers(SSc_Nor,only.pos = T ,min.pct = 0.1 ,logfc.threshold = 0.25)
top10 <- SSc_Nor.markers%>%group_by(cluster)%>%top_n(n=10, wt=avg_log2FC)
find <- top10[top10$gene == 'MS4A1', ]
find <- SSc_Nor.markers[SSc_Nor.markers$gene == "CD19", ]
DoHeatmap(SSc_Nor,features = top10$gene)
FeaturePlot(SSc_Nor,features = c("CD4"),reduction = 'umap',label = TRUE)
FeaturePlot(SSc_Nor,features = c("CD4"),reduction = 'umap')
FeaturePlot(SSc_Nor,features = c("XIST"),split.by = 'orig.ident',reduction = 'umap')
FeaturePlot(SSc_Nor,features = c("KRT5"),reduction = 'umap')
FeaturePlot(SSc_Nor,features = c("KRT5"),reduction = 'umap')
# 定义标记基因
marker_genes <- list(
  Fibroblast = c('FAP', 'S100A4', 'PDGFRB'),
  Pericyte = c('PDGFRB', 'ACTA2', 'RGS5'),
  Smooth_Muscle = c('ACTA2', 'MYH11', 'CNN1'),
  Keratinocyte = c('KRT14', 'KRT5', 'DSG1'),
  Nerve = c('TUBB3', 'MAP2', 'NEFM'),
  Mast_Cell = c('KIT', 'FCER1A', 'TPSAB1'),
  B_Cell = c( 'CD79A', 'CD19'),
  T_Cell = c( 'CD4', 'CD8A'),
  Melanocyte = c('MLANA', 'TYR', 'DCT'),
  Endothelial = c('PECAM1', 'VWF', 'CDH5'),
  Eccrine_Gland = c('KRT19', 'AQP5', 'KRT7'),
  Myeloid = c('ITGAM', 'CD14', 'CSF1R')
)
all_marker_genes <- unlist(marker_genes)
# 将默认assay设置为RNA
DefaultAssay(SSc_Nor) <- "RNA"

# 使用 DotPlot 进行可视化
DotPlot(SSc_Nor, features = B_T_cell, group.by = 'seurat_clusters', dot.scale = 6) +
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = round(avg.exp, 2)), size = 2)

DefaultAssay(SSc_Nor) = "RNA"
DotPlot(SSc_Nor,features = unique(marker_genes),group.by = 'RNA_snn_res.1') + coord_flip()
expr_data <- SSc_Nor@assays$RNA@layers$counts
allmarkers <- FindAllMarkers(expr_data,logfc.threshold = 0.5,min.pct = 0.1,only.pos = T)

SSc_Nor <- RunAzimuth(
  query = SSc_Nor, 
  reference = "skinref", 
  reference.assay = "RNA", 
  query.assay = "RNA", 
  normalization.method = "SCT", 
  dims = 1:30
)
DimPlot(SSc_Nor, reduction = "umap", group.by = c("orig.ident", "seurat_annotations"))
#---------------------------------------------------------------------------------------------------------
setwd('H:/bioinformatics/transcriptomics/SSc scRNA-RNAseq/SSc数据/GSE195452_scRNA_skin_blood/GSE195452_scRNA_skin_blood')
txt_files <- list.files(pattern = "*.txt", full.names = TRUE)
sample_names <- sub("\\./", "", sub("_.*", "", txt_files))
seurat_list <- list()
for (i in 88:length(txt_files)) {
  cat("Processing file", i, "of", length(txt_files), "\n")
  
  # 读取每个样本的表达矩阵
  counts <- fread(txt_files[i])
  
  # 使用第一列作为行名
  counts <- column_to_rownames(counts, var = colnames(counts)[1])
  
  # 替换行名中的下划线为破折号
  rownames(counts) <- gsub("_", "-", rownames(counts))
  counts_matrix <- as.matrix(counts)
  # 创建Seurat对象
  seurat_obj <- CreateSeuratObject(counts = as.matrix(counts_matrix),
                                   project = sample_names[i],
                                   min.features = 200, 
                                   min.cells = 3)
  
  # 计算线粒体基因比例
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # 过滤细胞
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  
  # 数据预处理
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  
  # 将Seurat对象添加到列表中
  seurat_list[[i]] <- seurat_obj
}



