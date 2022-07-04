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

prepare_mean_dna <- function(seurat_scDNA){
  
  ## Filter the noise cells and S phase cells
  if(length(which(seurat_scDNA$celltype == 'normal')) == 0){
    seurat_scDNA_subset <- seurat_scDNA
  }else{
    cells = union(which(seurat_scDNA$celltype == 'normal'), which(seurat_scDNA$celltype == 'tumor'))
    seurat_scDNA_subset <- subset(seurat_scDNA, cells = cells)
  }

 
  CNVmatrix <- seurat_scDNA_subset@assays$RNA@counts
  label <- as.character(seurat_scDNA_subset$subclone)
  
  if('normal' %in% label){
    label[which(label == 'normal')] <- 'Normal'
  }
  
  dna_df <- dplyr::as_data_frame(t(CNVmatrix)) %>%
    dplyr::mutate(label = label) %>% dplyr::group_by(label) %>%
    dplyr::summarise_all(mean) %>% dplyr::ungroup() 
  
  subclone <- dna_df$label
  
  dna_df <- dna_df %>% dplyr::select(-c(1)) 
  
  mean_dna_df <- dplyr::as_data_frame(t(dna_df)) %>%
    dplyr::mutate(chr = scdna_matrix_locs$chr) %>%
    dplyr::group_by(chr) %>% dplyr::summarise_all(mean) %>%
    dplyr::ungroup() %>% 
    dplyr::filter(chr != 'chrX' & chr !='chrY') %>%
    dplyr::arrange(as.numeric(gsub('chr', '', chr))) %>% 
    dplyr::select(-c(1)) %>% t()
  
  rownames(mean_dna_df) <- subclone
  colnames(mean_dna_df) <- paste0('chr', 1:22)  
  
  return(mean_dna_df)
}

prepare_mean_rna <- function(seurat_obj){
  
  label <- as.character(seurat_obj$corrlabel)
  if('normal' %in% label){
    label[which(label == 'normal')] <- 'Normal'
  }
  
  # normalized data
  mat = seurat_obj@assays$RNA@data
  common_gene = intersect(chr_genes$gene, rownames(mat))
  mat = mat[common_gene, ]
  
  chr_df <- chr_genes %>% filter(gene %in% common_gene) %>% filter(!duplicated(gene))
  
  ## Z-score for mat
  # mat_zscore <- apply(mat, 1, function(x){
  #   if(sd(x) == 0){
  #   s = 1}else{
  #   s = sd(x)
  #   }
  #   return((x-mean(x))/s)
  #   })
  
  ## Calculate the mean of expression on subclones
  rna_df <- dplyr::as_data_frame(t(mat)) %>%
    dplyr::mutate(label = label) %>% dplyr::group_by(label) %>%
    dplyr::summarise_all(mean) %>% dplyr::ungroup() 
  
  # ## Calculate the mean of expression on subclones
  # rna_df <- dplyr::as_data_frame(mat_zscore) %>%
  #   dplyr::mutate(label = label) %>% dplyr::group_by(label) %>%
  #   dplyr::summarise_all(mean) %>% dplyr::ungroup() 
  
  subclone <- rna_df$label
  
  rna_df <- rna_df %>% dplyr::select(-c(1)) 
  
  mean_rna_df <- dplyr::as_data_frame(t(rna_df)) %>%
    dplyr::mutate(chr = chr_df$chr) %>%
    dplyr::group_by(chr) %>% dplyr::summarise_all(mean) %>%
    dplyr::ungroup() %>% 
    dplyr::filter(chr != 'chrX' & chr !='chrY') %>%
    dplyr::arrange(as.numeric(gsub('chr', '', chr))) %>% 
    dplyr::select(-c(1)) %>% t()
  
  rownames(mean_rna_df) <- subclone
  colnames(mean_rna_df) <- paste0('chr', 1:22)  
  
  return(mean_rna_df)
}



