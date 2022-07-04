library(ggplot2)
library(optparse)
library(future)
library(optparse)
library(Seurat)

library(IRanges)
library(data.table)
library(dplyr)
library(ggplot2)
library(mixtools)
library(GenomicRanges)
library(RColorBrewer)

source('/mnt/ix1/Projects/M070_200622_GI_multiomics/scDNA/00_Code/code_v1/Plot_all_chr_heatmap_latest.R', echo=TRUE)
source('/mnt/ix1/Projects/M070_200622_GI_multiomics/scDNA/00_Code/code_v1/Filter_noise.R', echo=TRUE)

options(future.globals.maxSize = 24000 * 1024^2)
plan('multiprocess', workers=24)


# setwd('/mnt/ix1/Projects/M070_200622_GI_multiomics/scDNA/A12_210825_seurat/P5931_tumor')

option_list <- list(
  make_option(c('-s', '--sample'), action='store', default = NULL, type='character', help = 'Sample name')
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$sample)) {
  print_help(opt_parser)
  stop("Need sample name.")
}

### Run seurat cluster for noise cells
run_seurat <- function(mat){
  
  row.names(mat) <- paste0('segement', 1:dim(mat)[1])
  dims.reduce <- 50
  ncells <- dim(mat)[2]
  
  r = 0.2
  # if ( ncells <= 600 ){
  #   r <- 0.1
  # }else if(600 < ncells & ncells < 1200){
  #   r <- 0.2
  # }else if (1200 <= ncells & ncells  < 1800){
  #   r <- 0.3
  # }else if (1800 <= ncells & ncells < 2400){
  #   r <- 0.4
  # }else{
  #   r <- 0.5
  # }
  
  y <- tryCatch({
    Iden_signal_segements(mat)
  },warning= function(w){
    # index_segement <- Iden_signal_segements(mat)
    # print('warning')
    return(0)
  }, error = function(e){
    # index_segement <- seq(1, dim(mat)[1])
    # print('error')
    return(0)
  })
  
  if(y == 0){
    index_segement <- seq(1, dim(mat)[1])
  }else{
    index_segement <- Iden_signal_segements(mat)
  }
  
  write.csv(scdna_matrix_locs$chr[index_segement], file='signal_chr.csv')
  seurat_scDNA <- CreateSeuratObject(counts = as.matrix(mat), project='seurat-v3', min.cells = -Inf, min.features = -Inf)
  
  signal_segement <- rownames(mat)[index_segement]
  seurat_scDNA<- subset(seurat_scDNA, features = signal_segement)
  seurat_scDNA <- seurat_scDNA %>% ScaleData(do.center=F, do.scale = F) %>% RunPCA(verbose=T, features=rownames(seurat_scDNA), npcs = min(50, ncol(seurat_scDNA)-1)) %>% 
    RunUMAP(umap.method = 'uwot', n.neighbors = min(20, ncol(seurat_scDNA)), dims = 1:min(20, ncol(seurat_scDNA)-1), verbose=F) %>% FindNeighbors(reduction='umap', dims=1:2)%>% FindClusters(reduction = 'umap', resolution = r)
  
  Idents(seurat_scDNA) <- paste0('C', as.numeric(Idents(seurat_scDNA)))
  return(seurat_scDNA)
}

### plot heatmap
plot_heatmap <- function(scdna_matrix_merge, seurat_scDNA, title){
  
  labels <- Idents(seurat_scDNA)
  
  sum_row <- apply(scdna_matrix_merge, 1, mean)
  row_index <- which(sum_row != 0)
  scdna_matrix_merge <- scdna_matrix_merge[row_index, ]
  scdna_matrix_locs <- scdna_matrix_locs[row_index, ]
  
  scdna_matrix_merge$chr <- scdna_matrix_locs$chr
  scdna_matrix_merge$bin <- scdna_matrix_locs$bin
  
  scdna_matrix_merge <- scdna_matrix_merge %>% arrange(bin)
  scdna_matrix_locs <- scdna_matrix_locs %>% arrange(bin)
  
  scdna_matrix_merge$chr <- NULL
  scdna_matrix_merge$bin <- NULL
  
  scdna_matrix_locs$seg <- paste0('seg', 1:dim(scdna_matrix_locs)[1])
  scdna_matrix_locs[['Gene_ID']] = as.character(scdna_matrix_locs[["seg"]])
  scdna_matrix_locs[['Gene_Chr']] = as.character(scdna_matrix_locs[["chr"]])
  row.names(scdna_matrix_locs) = scdna_matrix_locs[['Gene_ID']]
  
  rownames(scdna_matrix_merge) <- scdna_matrix_locs$seg
  
  gene_chr = rle(scdna_matrix_locs[["Gene_Chr"]][match(rownames(scdna_matrix_merge), scdna_matrix_locs[["Gene_ID"]])])
  gene_chr_mark = c(0, cumsum(gene_chr$lengths))[seq_along(gene_chr$lengths)]+1
  names(gene_chr_mark) = gene_chr$values
  
  mat = scdna_matrix_merge
  label <- as.numeric(Idents(seurat_scDNA))
  
  # Recorder cells in each subclones based on the distance to normal CNV profile
  normal_CNV <- matrix(2, nrow = 1, ncol = dim(scdna_matrix_merge)[1])
  d <- as.matrix(dist(t(cbind(t(as.matrix(normal_CNV)), mat))))[1, ]
  d <- d[1:dim(mat)[2]+1]
  
  mat <- SortCells(mat, label, d)
  mat[mat > 6] = 6
  
  labels <- rep(unique(Idents(seurat_scDNA)), table(sort(label)))
  labels <- data.frame(cluster=labels)
  # Do not need to reorder the subclones
  # obj <- Sort_subclones(mat, labels, hc)
  # mat <- obj[[1]]
  # labels <- obj[[2]]
  Plot_CNV_heatmap(mat, labels, gene_chr_mark, title)
}


project_path <- '/mnt/ix1/Projects/M070_200622_GI_multiomics/scDNA'
load('scdna_matrix_50k_all_barcodes.Robj')
load('scdna_matrix_locs.Robj')
load('seurat_scDNA.Robj')

## cell metrics 
metrics_file <- paste0(project_path, '/A00_cellranger_dna/', opt$sample, '/outs/per_cell_summary_metrics.csv')
per_cell_metrics <- read.csv(metrics_file, sep=",", header=T)

# Identify the technical noise cells 
num_NA <- apply(scdna_matrix_merge_allbarcodes, 2, function(x){length(which(x==0))})
prop_NA <- num_NA/dim(scdna_matrix_merge_allbarcodes)[1]
prop_format <- data.frame(prop = prop_NA, index = seq(1, length(prop_NA)))

p<- ggplot(data=prop_format, aes(x= index, y=prop)) + 
  geom_point(size = 1, alpha = 0.5) + 
  xlab("Cell index") +
  ylab("CNV dropout proportion") +
  scale_y_continuous(breaks= seq(0, max(prop_format$prop), 0.1)) + 
  geom_hline(aes(yintercept = 0.1, colour = '#990000', linetype = 'dashed'), show.legend = FALSE) + 
  theme(plot.title = element_text(size = 11, color = 'black', face = 'bold', hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 15, color = 'black', face = 'bold')) +
  theme(axis.text.y = element_text(size = 15, color = 'black', face = 'bold')) +
  theme(axis.title.x = element_text(size = 15, color = 'black', face = 'bold')) + 
  theme(axis.title.y = element_text(size = 15, color = 'black', face = 'bold')) 
  # theme(legend.text= element_text(size=15,color="black", face= "bold", vjust=0.5, hjust=0.5)) +
  # theme(legend.title = element_text(size=15,color="black", face= "bold"))
ggsave(p, file = 'NA_proportion.pdf')

technoise_index <- which(prop_NA > 0.1)
label_technoise <- matrix(0, nrow = dim(seurat_scDNA)[2])
label_technoise[technoise_index] <- 1

write.csv(table(per_cell_metrics$is_noisy[technoise_index]), file = 'technoise.csv')

scdna_matrix_merge_technoise <- scdna_matrix_merge_allbarcodes[, technoise_index]
# save(scdna_matrix_merge_technoise, file='scdna_matrix_technoise.Robj')
seurat_scDNA$is.noisy <- per_cell_metrics$is_noisy
seurat_scDNA$is.technoisy <- label_technoise

celltype <- matrix(0, nrow=dim(seurat_scDNA)[2])
celltype[which(per_cell_metrics$is_noisy == 1)] <- 'cellranger noise'
celltype[which(label_technoise == 1)] <- 'technical noise'

# Filter the technical noisy cells
seurat_scDNA_sub1 <- subset(seurat_scDNA, cells = per_cell_metrics$barcode[which(label_technoise == 0)])
# Filter the cellranger-cnv noise
seurat_scDNA_sub2 <- subset(seurat_scDNA_sub1, cells = per_cell_metrics$barcode[which(per_cell_metrics$is_noisy == 0)])

### Identify the normal cells
tmp_matrix <- as.data.frame(seurat_scDNA_sub2@assays$RNA@counts)
tmp_matrix[tmp_matrix == 0] <- NA
tmp_matrix$chr <- scdna_matrix_locs$chr

cell_averages <- tmp_matrix %>% group_by(chr) %>% summarise_all(list(mean), na.rm=T) %>% as.data.frame()
cell_averages$chr <- NULL
cell_averages_melt <- data.table::melt(cell_averages, measure.vars=colnames(cell_averages))
detach(package:plyr)
cell_averages_classification <- cell_averages_melt %>% dplyr::group_by(variable) %>% summarise(c=any(value > 2.5, na.rm=T))

celltype[which(celltype == 0)] <- cell_averages_classification$c
celltype[which(celltype == "FALSE")] <- 'normal'

## Filter the normal cells
seurat_scDNA_tumor <- subset(seurat_scDNA, cells = colnames(seurat_scDNA_sub2)[cell_averages_classification$c])

## Identify the replication cells
locs <- scdna_matrix_locs %>% filter(chr != 'chrX' & chr != 'chrY')
mat <- seurat_scDNA_tumor@assays$RNA@counts[1:length(locs$chr), ]
# Recorder cells in each subclones based on the distance to normal CNV profile
normal_CNV <- matrix(2, nrow = 1, ncol = dim(mat)[1])
d <- as.matrix(dist(t(cbind(t(as.matrix(normal_CNV)), mat))))[1, ]
d <- d[1:dim(mat)[2]+1]
# Filter the noisy cells by EM (mixture normal distributions)
# start.par <- mean(d, na.rm = TRUE) + sd(d, na.rm = TRUE) * runif(2)
estimate_mix_model <- normalmixEM(d)
est <- data.frame(lambda = estimate_mix_model$lambda, mu = estimate_mix_model$mu, sigma = estimate_mix_model$sigma)
write.csv(est, file = 'est_mix_model.csv')

pdf('mixture_cells.pdf')
mixtools::plot.mixEM(estimate_mix_model, whichplots = 2, main2 = 'S phase cells', xlab2 = 'The Euclidean distance to normal cells' , lwd2 = 3, marginal = TRUE)
dev.off()

#S phase cells from two categories, either from the distribution with big mu and big sigma, or
#the distribution from small mu with big sigma
mu <- estimate_mix_model$mu
lambda <- estimate_mix_model$lambda
sigma <- estimate_mix_model$sigma
index_label <- apply(estimate_mix_model$posterior, 1, function(x){which(x == max(x))})
if(max(mu) > 2*min(mu)){
  index_Sphase <- which(mu == max(mu))
}else{
  index_Sphase <- which(sigma == max(sigma))
}

index_Scells <- which(index_label == index_Sphase)

index_label[index_Scells] <- 'S phase'
index_label[which(index_label != 'S phase')] <- 'tumor'
celltype[which(celltype == 'TRUE')] <- index_label
seurat_scDNA$celltype <- celltype

## Run Seurat for tumor cells
seurat_scDNA_subset <- subset(seurat_scDNA, cells = which(celltype == 'tumor'))
seurat_scDNA_subset <- run_seurat(seurat_scDNA_subset@assays$RNA@counts)

celltype <- seurat_scDNA$celltype
celltype[which(celltype == 'tumor')] <- paste0('C', as.numeric(seurat_scDNA_subset$seurat_clusters))
seurat_scDNA$subclone <- celltype

save(seurat_scDNA, file = 'seurat_scDNA.Robj')
save(seurat_scDNA_subset, file = 'seurat_scDNA_tumor.Robj')


# ################
# scdna_matrix_merge_technoise$chr <- scdna_matrix_locs$chr
# scdna_matrix_merge_technoise$bin <- scdna_matrix_locs$bin
# scdna_matrix_merge_technoise <- scdna_matrix_merge_technoise %>% filter(chr != 'chrX' & chr != 'chrY')
# scdna_matrix_locs <- scdna_matrix_locs %>% filter(chr != 'chrX' & chr != 'chrY')
#
# scdna_matrix_merge_technoise$chr <- NULL
# scdna_matrix_merge_technoise$bin <- NULL
#
# ### The heatmap of the cells with a lot of zero
# label_technoise <- as.factor(rep('noisy cells', dim(scdna_matrix_merge_technoise)[2]))
# seurat_technoise <- CreateSeuratObject(counts = scdna_matrix_merge_technoise, project='seurat-v3', min.cells = -Inf, min.features = -Inf)
# Idents(seurat_technoise) <- label_technoise
# plot_heatmap(scdna_matrix_merge_technoise, seurat_technoise, './tech_noise/delete_cells_heatmap.pdf')



