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

mean_dna_df <- readRDS('mean_dna_df.rds')
mean_rna_df <- readRDS('mean_rna_df.rds')
sig.chr <- read.table('par_chr.txt')
sig.chr$V1 <- gsub('C', 'c', sig.chr$V1)

mean_dna_df <- round(mean_dna_df[ , sig.chr$V1])
mean_rna_df <- mean_rna_df[ , sig.chr$V1]

if(is.null(rownames(mean_rna_df))){
  mean_rna_df <- as.matrix(t(mean_rna_df))
  rownames(mean_rna_df) <- 'C1'
}

if(dim(mean_rna_df)[1] == dim(mean_dna_df)[1]){
  mean_dna_df <- mean_dna_df[rownames(mean_rna_df), ]
  if(identical(rownames(mean_dna_df), rownames(mean_rna_df))){
    print('TRUE')
  }else{
    print('Reorder the subclones between scRNA and scDNA')
  }
}else{
  mean_dna_df <- mean_dna_df[rownames(mean_rna_df), ]
  identical(rownames(mean_dna_df), rownames(mean_rna_df))
}

if('Normal' %in% rownames(mean_dna_df)){
  mean_dna_df <- mean_dna_df[-grep("Normal",rownames(mean_dna_df)), ]
  if(is.null(dim(mean_dna_df))){
    mean_dna_df <- as.matrix(t(mean_dna_df))
    rownames(mean_dna_df) <- 'C1'
  }
}

if('Normal' %in% rownames(mean_rna_df)){
  mean_rna_df <- mean_rna_df[-grep("Normal",rownames(mean_rna_df)), ]
  if(is.null(dim(mean_rna_df))){
    mean_rna_df <- as.matrix(t(mean_rna_df))
    rownames(mean_rna_df) <- 'C1'
  }
}

if(is.null(rownames(mean_dna_df))){
  mean_dna_df <- as.matrix(t(mean_dna_df))
  rownames(mean_dna_df) <- 'C1'
}

if(is.null(rownames(mean_rna_df))){
  mean_rna_df <- as.matrix(t(mean_rna_df))
  rownames(mean_rna_df) <- 'C1'
}

### Start to plot
# f1 = circlize::colorRamp2(seq(min(mean_dna_df), max(mean_dna_df), length = 3),
#                           c('#3182bd', '#eff3ff',  '#e6550d'))

if(min(mean_dna_df) == 2){
  f1 = circlize::colorRamp2(seq(2, max(mean_dna_df), length = 4),
                            c('#eff3ff', '#fdbe85', '#fd8d3c', '#e6550d'))
}else{
  f1 = circlize::colorRamp2(c(min(mean_dna_df), 2, max(mean_dna_df)),
                            c('#3182bd', '#eff3ff', '#e6550d'))
}

f2 = circlize::colorRamp2(seq(min(mean_rna_df), max(mean_rna_df), length = 3),
                          c("#67A9CF", "#F7F7F7", "#EF8A62"))

htdna <- ComplexHeatmap::Heatmap(mean_dna_df,
                                 col = f1,
                                 name = 'avg_CNV',
                                 # border = TRUE,
                                 rect_gp = gpar(col='white', lwd=1),
                                 column_title = 'Average of Copy number',
                                 cluster_columns = FALSE,
                                 cluster_rows = FALSE)

htrna <- ComplexHeatmap::Heatmap(mean_rna_df,
                                 col = f2,
                                 name = 'avg_expression',
                                 # border = TRUE,
                                 rect_gp = gpar(col='white', lwd=1),
                                 column_title = 'Average of gene expression',
                                 cluster_columns = FALSE,
                                 cluster_rows = FALSE)
  
### If there is only one subclone
if(identical(rownames(mean_dna_df), rownames(mean_rna_df)) && dim(mean_dna_df)==1){
  ht = htrna %v% htdna
}else if(identical(rownames(mean_dna_df), rownames(mean_rna_df)) && dim(mean_dna_df) > 1){
  ht = htrna + htdna
}

pdf('Combine_average_subclones.pdf', width = 8, height = 4)
ComplexHeatmap::draw(ht,
                     ht_gap = unit(0.4, 'cm'),
                     heatmap_legend_side = 'left')  
dev.off()
  

  





