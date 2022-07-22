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
library(ggsignif)

mean_dna_df <- readRDS('mean_dna_df.rds')
mean_rna_df <- readRDS('mean_rna_df.rds')

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

dna_long <- melt(mean_dna_df)
names(dna_long) <- c('subclone', 'chr', 'Avg_CNV')
dna_long$Avg_CNV <- round(dna_long$Avg_CNV)

rna_long <- melt(mean_rna_df)
names(rna_long) <- c('subclone', 'chr', 'Avg_expression')

remin_subclone <- intersect(unique(rna_long$subclone), unique(dna_long$subclone))
dna_long <- dna_long %>% filter(subclone %in% remin_subclone) 
dna_long$subclone <- as.character(dna_long$subclone)
rna_long <- rna_long %>% filter(subclone %in% remin_subclone)
rna_long$subclone <- as.character(rna_long$subclone)

identical(dna_long$subclone, rna_long$subclone)
identical(rna_long$chr, rna_long$chr)

final_long <- dna_long %>% dplyr::mutate(Avg_exp = rna_long$Avg_expression)
## Filter the CNV which one have one points

tmp_cnv <- plyr::count(final_long$Avg_CNV)$x[which(plyr::count(final_long$Avg_CNV)$freq > 1)]

final_long <- final_long %>% filter(Avg_CNV %in% tmp_cnv)

plot_violin_cnv <- function(final_long){
  
  p <- ggplot(final_long, aes(x = Avg_CNV, y = Avg_exp, fill = as.factor(Avg_CNV))) + 
    geom_violin(trim = FALSE) +
    geom_signif(comparsion = list(c('1', '2'), c('2', '3')), map_signif_level = TRUE) + 
    # geom_jitter(shape = 16, position = position_jitter(0.1))
    scale_fill_brewer(palette="Set2") + 
    stat_summary(fun.data="mean_sdl", mult=0.5, width = 0.1,
                 geom="crossbar", color = 'black') + 
    theme(
      axis.title.x = element_text(size = 15, color = 'black', face = 'bold', hjust = 0.5),
      axis.title.y = element_text(size = 15, color = 'black', face = 'bold', hjust = 0.5),
      axis.text.y = element_text(size = 12, color = 'black', face = 'bold', hjust = 0.5),
      strip.text = element_text(size = 15, color = 'black', face = 'bold', hjust = 0.5),
      legend.title = element_text(size = 15, color = 'black', face = 'bold', hjust = 0.5),
      legend.text = element_text(size = 10, color = 'black', face = 'bold', hjust = 0.5),
      axis.text.x = element_text(size = 12, angle = 45, color = 'black', face = 'bold', hjust = 1),
      axis.line.x = element_line(size = 0.5, colour = 'black'), #加上x轴的坐标轴
      axis.line.y = element_line(size = 0.5, colour = 'black'),
      panel.background = element_blank(),
      legend.position="none"
    ) + 
    ylab('Average chromosomal expresssion of subclones') +
    xlab('Average copy number of subclones') 
  # guides(colour = guide_legend(title = 'Subclone'))
  
  ggsave(p, file = 'violin_plot_cnv.pdf', width = 8, height = 6)
}

plot_violin_cnv(final_long)



