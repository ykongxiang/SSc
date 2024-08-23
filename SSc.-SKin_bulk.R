library(data.table)
library(dplyr)     
library(ggplot2)
setwd('H:/bioinformatics/transcriptomics/SSc scRNA-RNAseq/SSc数据/GSE254716_gene_expected_count.annot.txt')
SSc_skin <- fread("H:/bioinformatics/transcriptomics/SSc scRNA-RNAseq/SSc数据/GSE254716_gene_expected_count.annot.txt/GSE254716_gene_expected_count.annot.txt")
SSc_skin <- SSc_skin[,c(-6,-7,-9,-11,-13)]
SSc_