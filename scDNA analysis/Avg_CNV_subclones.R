library(SummarizedExperiment)
library(Seurat)
library(dplyr)
library(annotables)
library(optparse)
library(future)

source('/home/xiangqi/shared/R_scripts/InitializeCNV.R')

options(future.globals.maxSize = 6000 * 1024^2)
plan('multiprocess', workers= 24)

load('seurat_scDNA.Robj')
CNVmatrix <- readRDS('CNV_gene_matrix.rds')
identical(colnames(CNVmatrix), colnames(seurat_scDNA))

cell_nonoise <- which(seurat_scDNA$celltype %in% c('normal', 'tumor'))
## Filter the noise cells
seurat_scDNA_subset <- subset(seurat_scDNA, cells = cell_nonoise)
CNVmatrix_subset <- as.data.frame(t(CNVmatrix[, cell_nonoise]))
CNVmatrix_subset$subclone <- seurat_scDNA_subset$subclone

CNV_subset_mean <- CNVmatrix_subset %>% dplyr::group_by(subclone) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) 

subclone <- CNV_subset_mean$subclone
CNV_subset_mean$subclone <- NULL
CNV_subset_mean <- round(CNV_subset_mean)

colnames(CNV_subset_mean) <- rownames(CNVmatrix)
rownames(CNV_subset_mean) <- subclone
  
saveRDS(CNV_subset_mean, file = 'CNV_subclones.rds')

