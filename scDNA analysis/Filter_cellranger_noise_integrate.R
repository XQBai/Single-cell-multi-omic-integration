library(pheatmap)
library(ggplot2)
library(Seurat)
library(stringr)
library(Rtsne)
#library(uwot)
library(mclust)
library(optparse)
library(biomaRt)
library(cowplot)
library(ggplotify)
library(RColorBrewer)
library(minerva)
library(data.table)
library(biomaRt)
library(future)

option_list <- list(
  # make_option(c("-i", "--input"), action="store", default=NULL, type='character', help="Input base folder for mtx and tsv files"),
  make_option(c("-s", "--sample"), action="store", default=NULL, type='character', help="Sample name/output prefix")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


if (is.null(opt$sample)) {
  print_help(opt_parser)
  stop("Need sample name.")
}

#removing the noisy cells based on cellranger's identify
project_path <- '/mnt/ix1/Projects/M070_200622_GI_multiomics/scDNA/A13_210910_mainfigure1/A04_210910_TCGA/A01_210910_armMatrix/'
seurat_path <- '/mnt/ix1/Projects/M070_200622_GI_multiomics/scDNA/A12_210825_seurat/'

## The identified noise cells and S phase cells
load(paste0(seurat_path, opt$sample, '/seurat_scDNA.Robj'))

## Arm matrix
Armfile <- paste0(project_path, opt$sample, '_CNV_Arm_matrix.Robj')
load(Armfile)
DNAmatrixArm <- DNAmatrixArm[, which(seurat_scDNA$cellcomponment != 'Noise')]
save(DNAmatrixArm, file=paste0(project_path, opt$sample,'_CNV_Arm_matrix_nonoise.Robj'))

# #### Gene Matrix
# Genefile <- paste0(project_path, opt$sample, '_CNV_gene_matrix.Robj')
# load(Genefile)
# CNV_gene_Matrix <- CNV_gene_Matrix[, which(per_cell_metrics$is_noisy == 0)]
# save(CNV_gene_Matrix, file = paste0(project_path, opt$sample,'_CNV_gene_matrix_nonoise.Robj'))








