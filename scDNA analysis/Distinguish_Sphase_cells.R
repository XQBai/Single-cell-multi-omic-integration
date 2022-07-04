library(Seurat)

source('/mnt/ix1/Projects/M070_200622_GI_multiomics/scDNA/00_Code/code_v1/Plot_all_chr_heatmap_latest.R', echo=TRUE)

load('scdna_matrix_locs.Robj')
load('seurat_scDNA.Robj')

### plot heatmap
plot_heatmap <- function(scdna_matrix_merge, seurat_scDNA, title){
  
  labels <- Idents(seurat_scDNA)
  
  scdna_matrix_merge$chr <- scdna_matrix_locs$chr
  
  scdna_matrix_merge <- scdna_matrix_merge %>% filter(chr != 'chrX' & chr != 'chrY') 
  scdna_matrix_merge$chr <- NULL
  scdna_matrix_locs <- scdna_matrix_locs %>% filter(chr != 'chrX' & chr != 'chrY') 
  
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
  Plot_CNV_heatmap(mat, labels, gene_chr_mark, title)
}

# cellranger_noise_mat <- seurat_scDNA@assays$RNA@counts[, which(seurat_scDNA$celltype == 'cellranger noise')]
# cnvmean <- apply(cellranger_noise_mat, 2, mean)
# celltype <- seurat_scDNA$celltype
# cellcycle <- celltype
# cellcycle[which(celltype == 'cellranger noise')][which(cnvmean < 2.5)] <- 'noise'
# cellcycle[which(celltype == 'cellranger noise')][which(cnvmean >= 2.5)] <- 'cellranger S phase'
# seurat_scDNA$cellcycle <- cellcycle
# index <- union(which(seurat_scDNA$celltype == 'cellranger noise'), which(seurat_scDNA$celltype == 'S phase'))
# seurat.noise <- subset(seurat_scDNA, cells = index)
# Idents(seurat.noise) <- seurat.noise$cellcycle
# plot_heatmap(as.data.frame(seurat.noise@assays$RNA@counts), seurat.noise, 'Subclones_split_S_phase.pdf')

### Distinguish the noise cells and S phase cells from the union of cellranger noise and S phase
index <- union(which(seurat_scDNA$celltype == 'cellranger noise'), which(seurat_scDNA$celltype == 'S phase'))
seurat.noise <- subset(seurat_scDNA, cells = index)
noise.mat <- seurat.noise@assays$RNA@counts
noise.mean <- apply(noise.mat, 2, mean)

identifylabel <- matrix(0, nrow = length(noise.mean), ncol = 1)
identifylabel[which(noise.mean >= 2.5)] <- 'S phase'
identifylabel[which(noise.mean < 2.5)] <- 'cellranger noise'
Idents(seurat.noise) <- identifylabel
plot_heatmap(as.data.frame(seurat.noise@assays$RNA@counts), seurat.noise, 'Subclones_split_noise_Sphase.pdf')

## Plot the heatmap of noise cells
cellcycle <- seurat_scDNA$celltype
cellcycle[index][which(noise.mean >= 2.5)] <- 'S phase'
cellcycle[index][which(noise.mean < 2.5)] <- 'cellranger noise'
seurat_scDNA$cellcycle <- cellcycle
cellcycle[which(cellcycle == 'technical noise')] = 'Noise'
cellcycle[which(cellcycle == 'cellranger noise')] = 'Noise'
cellcycle[which(cellcycle == 'tumor')] = 'G0/G1 phase'
cellcycle[which(cellcycle == 'normal')] = 'Normal'
seurat_scDNA$cellcomponment <- cellcycle
save(seurat_scDNA, file = 'seurat_scDNA.Robj')


index_noise <- union(which(seurat_scDNA$cellcycle == 'cellranger noise'), which(seurat_scDNA$cellcycle == 'technical noise'))
seurat.noise <- subset(seurat_scDNA, cells = index_noise)
noise.mat <- seurat.noise@assays$RNA@counts
Idents(seurat.noise) <- seurat.noise$cellcycle
plot_heatmap(as.data.frame(seurat.noise@assays$RNA@counts), seurat.noise, 'Subclones_technical_noise.pdf')

## Plot the heatmap of normal cells and S phase cells 
seurat.sphase <- subset(seurat_scDNA, cells = which(seurat_scDNA$cellcycle == 'S phase'))
sphase.mat <- seurat.sphase@assays$RNA@counts
Idents(seurat.sphase) <- seurat.sphase$cellcycle
plot_heatmap(as.data.frame(seurat.sphase@assays$RNA@counts), seurat.sphase, 'Subclones_S_phase.pdf')



