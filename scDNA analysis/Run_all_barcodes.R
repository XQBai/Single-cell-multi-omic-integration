library(Seurat)
library(data.table)
#library(dbscan)
library(dplyr)
library(Matrix)
library(mixtools)
library(optparse)
library(future)

source('/mnt/ix1/Projects/M070_200622_GI_multiomics/scDNA/00_Code/Filter_noise.R', echo=TRUE)

options(future.globals.maxSize = 24000 * 1024^2)
plan('multiprocess', workers=24)

option_list <- list(
  make_option(c('-s', '--sample'), action='store', default = NULL, type='character', help = 'Sample name')
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$sample)) {
  print_help(opt_parser)
  stop("Need sample name.")
}

print(opt)
print(getwd())
sample <- opt$sample
#setwd('/mnt/ix1/Projects/M070_200622_GI_multiomics/scDNA/A03_200901_seurat/P5931_tumor')
#sample <- 'P5931_tumor'

project_path <- '/mnt/ix1/Projects/M070_200622_GI_multiomics/scDNA'
genome_reference <- fread(file.path(project_path, 'A02_bed_intersect', 'GRCh38_cellranger_20k.canonical.rownumbers.bed'), header=F, data.table=F)

CNV_file <- file.path(project_path, 'A00_200901_CNV', sample)
scdna_matrix <- fread(Sys.glob(file.path(CNV_file, '*.tsv'))[1])
barcodes <- read.csv(Sys.glob(file.path(CNV_file, '*.txt')), sep='\n', header = FALSE)
colnames(scdna_matrix) <- as.character(barcodes$V1)

bin_size <- 50 #number of 20kb bins

genome_reference <- genome_reference %>% group_by(V1) %>% mutate(bin=floor(row_number(V1)/bin_size))
genome_reference$bin_corrected <- cumsum(c(0,as.numeric(diff(genome_reference$bin))!=0))

scdna_matrix$chr <- genome_reference$V1
scdna_matrix$bin <- genome_reference$bin_corrected

#coarse graining
scdna_matrix_merge <- scdna_matrix %>% group_by(chr,bin) %>% summarise_all(list(median))
scdna_matrix_locs <- scdna_matrix_merge[,c("chr","bin")]
scdna_matrix_merge$chr <- NULL
scdna_matrix_merge$bin <- NULL

row.names(scdna_matrix_merge) <- paste0('segement', 1:dim(scdna_matrix_merge)[1])
dims.reduce <- 50
seurat_scDNA <- CreateSeuratObject(counts = as.matrix(scdna_matrix_merge), project='seurat-v3', min.cells = -Inf, min.features = -Inf)

seurat_scDNA <- seurat_scDNA %>% ScaleData(do.center=F, do.scale = F) %>% RunPCA(verbose=T, features=rownames(seurat_scDNA)) %>% RunUMAP(umap.method = 'uwot',n.neighbors = min(20, ncol(seurat_scDNA)), dims=1:dims.reduce) %>% 
  FindNeighbors(reduction='umap', dims=1:2)

scdna_matrix_merge_allbarcodes <- scdna_matrix_merge
save(scdna_matrix_merge_allbarcodes, file='scdna_matrix_50k_all_barcodes.Robj')
save(scdna_matrix_locs, file='scdna_matrix_locs.Robj')
save(seurat_scDNA, file='seurat_scDNA.Robj')






