library(dplyr)
library(GenomicRanges)
library(Seurat)
library(ggplot2)
library(future)

expr.infercnv.17 <- read.table('./test/expr.infercnv.17_HMM_predHMMi6.hmm_mode-samples.dat')
pred_cnv_genes <- read.table('./test/17_HMM_predHMMi6.hmm_mode-samples.genes_used.dat', header = T)
pred_cnv_regions <- read.table('./test/17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_regions.dat', header = T)
seurat <- readRDS('seurat_epi.rds')

label <- seurat$corrlabel
cellbarcode <- gsub('-', '.', colnames(expr.infercnv.17))
identical(cellbarcode, colnames(expr.infercnv.17))

###############
scdna_matrix_merge <- expr.infercnv.17
scdna_matrix_locs <- pred_cnv_genes
scdna_matrix_locs$chr <- paste0('chr', scdna_matrix_locs$chr)
scdna_matrix_merge$chr <- scdna_matrix_locs$chr
chrs <- paste0('chr', 1:22)

scdna_matrix_merge$chr <- NULL

sum_row <- apply(scdna_matrix_merge, 1, mean)
row_index <- which(sum_row != 0)
scdna_matrix_merge <- scdna_matrix_merge[row_index, ]
scdna_matrix_locs <- scdna_matrix_locs[row_index, ]

## Calculate the mean of expression on subclones
rna_df <- dplyr::as_data_frame(t(scdna_matrix_merge)) %>%
  dplyr::mutate(label = label) %>% dplyr::group_by(label) %>%
  dplyr::summarise_all(mean) %>% dplyr::ungroup() 

subclone <- rna_df$label
rna_df <- rna_df %>% dplyr::select(-c(1)) 

mean_rna_df <- dplyr::as_data_frame(t(rna_df)) %>%
  dplyr::mutate(chr = scdna_matrix_locs$chr) %>%
  dplyr::group_by(chr) %>% dplyr::summarise_all(mean) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(chr %in% paste0('chr', 1:22)) %>%
  dplyr::arrange(as.numeric(gsub('chr', '', chr))) %>% 
  dplyr::select(-c(1)) %>% t()

mean_rna_df <- round(mean_rna_df) - 1
rownames(mean_rna_df) <- subclone
colnames(mean_rna_df) <- paste0('chr', 1:22)

saveRDS(mean_rna_df, file = 'mean_rna_infercnv.rds')




