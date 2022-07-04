library(Seurat)
library(dplyr)
# library(bio3d)

expr.infercnv.17 <- read.table('./test/expr.infercnv.17_HMM_predHMMi6.hmm_mode-samples.dat')
pred_cnv_genes <- read.table('./test/17_HMM_predHMMi6.hmm_mode-samples.genes_used.dat', header = T)
pred_cnv_regions <- read.table('./test/17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_regions.dat', header = T)
seurat <- readRDS('seurat_epi.rds')

scdna_matrix_merge <- expr.infercnv.17 - (median(as.matrix(expr.infercnv.17)) - 2)
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
rna_df <- dplyr::as_data_frame(scdna_matrix_merge) %>%
  dplyr::mutate(chr = scdna_matrix_locs$chr) %>%
  dplyr::group_by(chr) %>% dplyr::summarise_all(mean) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(chr %in% paste0('chr', 1:22)) %>%
  dplyr::arrange(as.numeric(gsub('chr', '', chr))) %>% 
  dplyr::select(-c(1)) %>% t()

rna_df <- round(rna_df) - 1
# rownames(mean_rna_df) <- subclone
colnames(rna_df) <- paste0('chr', 1:22)

saveRDS(rna_df, file = 'rna_infercnv_mean_chrom.rds')

dna <- readRDS('mean_dna_df.rds')
dna <- round(dna)
label <- as.character(seurat$corrlabel)
label[which(label == 'normal')] <- 'Normal'
rna <- t(rna_df)

rsme <- function(rna, label, dna){
  
  rmsd <- matrix(0, nrow = length(unique(label)))
  for (i in 1:length(unique(label))){
    
    s <- as.character(unique(label)[i])
    index <- which(label == s)
    sc_rna <- rna[, index]
    sub_dna <- dna[s, ]
    
    ones <- matrix(1, nrow = 1, ncol = dim(sc_rna)[2])
    diff_mat <- as.matrix(sc_rna - as.matrix(sub_dna) %*% ones)
    l2_norm <- apply(diff_mat, 2, function(x){norm(x, type = '2')})
    rmsd[i] <- sqrt(sum(l2_norm)/length(index))
    # rmsd[i] <- sqrt(sum(l2_norm))
  }
  rmsd_avg <- mean(rmsd)
  
  return(rmsd_avg)
}

sig_chr <- c('chr3', 'chr7', 'chr8', 'chr20', 'chr21')
rna_signal <- rna[sig_chr, ]
dna_signal <- dna[, sig_chr]
rmsd_avg <- rsme(rna_signal, label, dna_signal)

write.csv(rmsd_avg, file = 'rmsd.csv')
