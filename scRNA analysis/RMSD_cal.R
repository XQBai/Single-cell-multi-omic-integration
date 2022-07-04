# Root-mean-square deviation
library(Seurat)
library(dplyr)
# library(bio3d)

expr.infercnv.17 <- read.table('./test/expr.infercnv.17_HMM_predHMMi6.hmm_mode-samples.dat')
pred_cnv_genes <- read.table('./test/17_HMM_predHMMi6.hmm_mode-samples.genes_used.dat', header = T)
pred_cnv_regions <- read.table('./test/17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_regions.dat', header = T)
seurat <- readRDS('seurat_epi.rds')

chr_genes <- read.table('/mnt/ix1/Projects/M070_200622_GI_multiomics/scRNA/00_Code/GSVA_tables/Chromo_pathway.csv', header=T, sep=",", stringsAsFactors=F)
sig.chr <- read.table('par_chr.txt')
## Find the common genes with signal genes
signal_genes <- chr_genes %>% filter(Pathway %in% sig.chr$V1)

label <- seurat$corrlabel
cellbarcode <- gsub('-', '.', colnames(expr.infercnv.17))
identical(cellbarcode, colnames(expr.infercnv.17))

CNV_subclone <- readRDS('CNV_subclones.rds')

common_gene <- intersect(rownames(expr.infercnv.17), colnames(CNV_subclone))
# common_gene <- intersect(common_gene, signal_genes$Genes)

rna <- expr.infercnv.17[common_gene, ] - (median(as.matrix(expr.infercnv.17))-2)
dna <- CNV_subclone[, common_gene]
rownames(dna) <- rownames(CNV_subclone)

identical(rownames(rna), colnames(dna))

# res <- rmsd(a, as.matrix(rna), ncore = 12)

rsme <- function(rna, label, dna){
  
  rmsd <- matrix(0, nrow = length(unique(label)))
  for (i in 1:length(unique(label))){
    
    s <- as.character(unique(label)[i])
    index <- which(label == s)
    sc_rna <- rna[, index]
    sub_dna <- dna[s, ]
    
    ones <- matrix(1, nrow = 1, ncol = dim(rna)[2])
    diff_mat <- as.matrix(sc_rna - t(sub_dna) %*% ones)
    l2_norm <- apply(diff_mat, 2, function(x){norm(x, type = '2')})
    rmsd[i] <- sqrt(sum(l2_norm)/length(index)/dim(sc_rna)[1])
    # rmsd[i] <- sqrt(sum(l2_norm))
  }
  rmsd_avg <- mean(rmsd)
  
  return(rmsd_avg)
}

rmsd_measure  <- rsme(rna, label, dna)
write.csv(rmsd_measure, file = 'rmsd.csv')

