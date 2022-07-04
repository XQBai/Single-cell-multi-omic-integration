library(Seurat)
library(optparse)
library(dplyr)
library(viridis)
library(cowplot)
library(future)

seurat_obj <- readRDS('seurat_epi.rds')
# label <- as.character(seurat_obj$corrlabel)
# 
# # # Cell annotation (use cluster# from original sample)
# cell_annot <- as.data.frame(label)
# cell_annot$barcode <- colnames(seurat_obj)
# write.table(cell_annot[,c(2,1)], 'idents.txt', sep="\t", row.names=F, col.names=F)

label <- as.character(seurat_obj$condition)

# # Cell annotation (use cluster# from original sample)
cell_annot <- as.data.frame(label)
cell_annot$barcode <- colnames(seurat_obj)
write.table(cell_annot[,c(2,1)], 'conditions.txt', sep="\t", row.names=F, col.names=F)


