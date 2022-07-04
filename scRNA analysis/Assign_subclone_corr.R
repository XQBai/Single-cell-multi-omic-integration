## Using corr to assign scRNA single cells into subclones 
library(SingleCellExperiment)
library(biomaRt)
library(Seurat)
library(dplyr)
library(Seurat)
library(future)
library(optparse)

options(future.globals.maxSize = 6000*1024^2)
plan('multiprocess', workers = 12)

# option_list = list(
#   make_option(c("-s", "--signal_chr"), action="store", default=NA, type='character',
#               help="signal chromosomes")
# )
# 
# opt = parse_args(OptionParser(option_list=option_list))

## Using the genes on chr3, chr7, chr8 and chr21 as input to clonealign 
chr_genes <- read.table('/mnt/ix1/Projects/M070_200622_GI_multiomics/scRNA/00_Code/GSVA_tables/Chromo_pathway.csv', header=T, sep=",", stringsAsFactors=F)
sig.chr <- read.table('par_chr.txt')

## Prepare the matrix
CNV_matrix <- readRDS('CNV_subclones.rds')
seurat_obj <- readRDS('seurat_epi.rds')
if (length(which(seurat_obj$Phase != 'S')) == 0){
  seurat_obj <- seurat_obj
}else{
  seurat_obj <- subset(seurat_obj, cells = which(seurat_obj$Phase != 'S'))
}

RNAscale <- seurat_obj@assays$RNA@data
RNAscale_tumor <- RNAscale[, which(seurat_obj$condition == 'tumor')]

## Find the common genes with signal genes
signal_genes <- chr_genes %>% filter(Pathway %in% sig.chr$V1)
common_genes <- intersect(intersect(rownames(RNAscale_tumor),
                                    colnames(CNV_matrix)), signal_genes$Genes)
subclone = rownames(CNV_matrix)
CNV_matrix <- CNV_matrix[which(subclone!= 'normal'), common_genes]
rownames(CNV_matrix) <- subclone[which(subclone != 'normal')]

RNAscale_tumor <- RNAscale_tumor[common_genes, ]

identical(colnames(CNV_matrix), rownames(RNAscale_tumor))

k = dim(CNV_matrix)[1]
corr_matrix <- matrix(0, nrow=dim(RNAscale_tumor)[2], ncol = k)

for (i in 1:k){
  print(i)
  mean_CNV <- t(CNV_matrix[i, ])
  corr_matrix[, i] <- as.matrix(apply(RNAscale_tumor, 2, function(x){cor(x, mean_CNV, method = 'pearson')[1]}))
}

label_corr <- apply(corr_matrix, 1, function(x){rownames(CNV_matrix)[which(x == max(x))[1]]})

#new label
label_corr_new <- seurat_obj$condition
label_corr_new[which(label_corr_new == 'tumor')] <- label_corr

seurat_obj$corrlabel <- factor(label_corr_new, levels = c('normal', rownames(CNV_matrix)))

saveRDS(seurat_obj, file = 'seurat_epi.rds')





