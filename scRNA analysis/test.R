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

if(rownames(mean_rna_df)=='C1' && 'normal' %in% rownames(mean_dna_df)){
  mean_dna_df <- as.matrix(t(mean_dna_df['C1', ]))
  rownames(mean_dna_df) <- 'C1'
}else if(dim(mean_dna_df)[1] == 2 && dim(mean_rna_df)[1] == 2){
  mean_dna_df <- t(as.matrix(mean_dna_df['C1', ]))
  mean_rna_df <- t(as.matrix(mean_rna_df['C1', ]))
  rownames(mean_dna_df) <- 'C1'
  rownames(mean_rna_df) <- 'C1'
}

# mean_dna_df <- mean_dna_df[rownames(mean_rna_df), ]
# mean_dna_df <- mean_dna_df[ , sig.chr]
# mean_rna_df <- mean_rna_df[ , sig.chr]

f1 = circlize::colorRamp2(seq(min(mean_dna_df), max(mean_dna_df), length = 3),
                          c('#3182bd', '#eff3ff',  '#e6550d'))
f2 = circlize::colorRamp2(seq(min(mean_rna_df), max(mean_rna_df), length = 3),
                          c("#67A9CF", "#F7F7F7", "#EF8A62"))

htdna <- ComplexHeatmap::Heatmap(mean_dna_df,
                                 col = f1,
                                 name = 'avg_CNV',
                                 # border = TRUE,
                                 rect_gp = gpar(col='white', lwd=1),
                                 row_title = 'scDNA',
                                 cluster_columns = FALSE,
                                 cluster_rows = FALSE)

htrna <- ComplexHeatmap::Heatmap(mean_rna_df,
                                 col = f2,
                                 name = 'avg_expression',
                                 # border = TRUE,
                                 rect_gp = gpar(col='white', lwd=1),
                                 row_title = 'scRNA',
                                 cluster_columns = FALSE,
                                 cluster_rows = FALSE)

ht = htrna %v% htdna

pdf('Combine_average_subclones.pdf', width = 6, height = 2)
ComplexHeatmap::draw(ht,
                     ht_gap = unit(0.4, 'cm'),
                     heatmap_legend_side = 'right')
dev.off()




