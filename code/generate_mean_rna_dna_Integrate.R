library(SingleCellExperiment)
library(data.table)
library(biomaRt)
library(dplyr)
library(Seurat)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(ggridges)
library(RColorBrewer)
library(ComplexHeatmap)
library(future)

source("/mnt/ix1/Projects/M070_200622_GI_multiomics/scRNA/00_Code/xiangqi/code_v1/prepare_matrix.R", echo=TRUE)

options(future.globals.maxSize = 6000 * 1024^2)
plan('multicore', workers= 18)

load('seurat_scDNA.Robj')
load('../scdna_matrix_locs.Robj')

seurat_obj <- readRDS('seurat_epi.rds')
chr_genes <- read.table('/mnt/ix1/Projects/M070_200622_GI_multiomics/scRNA/00_Code/xiangqi/bedfile/gene_locs_seg.tsv', 
                        header=F, sep="\t", stringsAsFactors=F)
names(chr_genes) <- c('gene', 'chr', 'start', 'end', 'seg') 

### DE genes between each subclones with normal clusters
Idents(seurat_obj) <- seurat_obj$corrlabel
label <- Idents(seurat_obj)

mean_dna_df <- prepare_mean_dna(seurat_scDNA)
mean_rna_df <- prepare_mean_rna(seurat_obj)
saveRDS(mean_dna_df, file = 'mean_dna_df.rds')
saveRDS(mean_rna_df, file = 'mean_rna_df.rds')

# sig.chr <- c('chr3', 'chr7', 'chr8', 'chr20', 'chr21')
# mean_dna_df <- mean_dna_df[ , sig.chr]
# mean_rna_df <- mean_rna_df[ , sig.chr]

# if('normal' %in% rownames(mean_dna_df)){
#   tmp_name <- setdiff(rownames(mean_dna_df), 'normal')
#   sort_name <- c('normal', tmp_name)
#   mean_dna_df <- mean_dna_df[sort_name,]
#   rownames(mean_dna_df) <- c('Normal', tmp_name)
#   rownames(mean_rna_df) <- c('Normal', tmp_name)
# }
# 
# f1 = circlize::colorRamp2(seq(min(mean_dna_df), max(mean_dna_df), length = 3),
#                           c('#3182bd', '#eff3ff',  '#e6550d'))
# f2 = circlize::colorRamp2(seq(min(mean_rna_df), max(mean_rna_df), length = 3),
#                           c("#67A9CF", "#F7F7F7", "#EF8A62"))
# 
# htdna <- ComplexHeatmap::Heatmap(mean_dna_df,
#                                  col = f1,
#                                  name = 'avg_CNV',
#                                  # border = TRUE,
#                                  rect_gp = gpar(col='white', lwd=1),
#                                  column_title = 'Average of Copy number',
#                                  cluster_columns = FALSE,
#                                  cluster_rows = FALSE)
# 
# htrna <- ComplexHeatmap::Heatmap(mean_rna_df,
#                                  col = f2, 
#                                  name = 'avg_expression',
#                                  # border = TRUE,
#                                  rect_gp = gpar(col='white', lwd=1),
#                                  column_title = 'Average of gene expression',
#                                  cluster_columns = FALSE,
#                                  cluster_rows = FALSE)
# 
# ht = htdna + htrna
# 
# pdf('Combine_average_subclones.pdf', width = 8, height = 4)
# ComplexHeatmap::draw(ht,
#                      ht_gap = unit(0.4, 'cm'),
#                      heatmap_legend_side = 'left')
# dev.off()





