library(Seurat)

pdf('umap_seurat_cluster.pdf')
DimPlot(seurat_obj, reduction = 'umap', label=TRUE)
dev.off()

## manually assign cell types to clusters
fib <- rep('Fibroblasts', 2)
end <- rep('Endothelial', 2)
macro <- c('Macrophages')
epi <- rep('Epithelial', 7)
den <- c('Dendritic')
Bpla <- rep('B plasma', 2)
tcell <- rep('T cells', 7)

celltypes <- c(fib, end, macro, epi, den, Bpla, tcell)

names(celltypes) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, celltypes)

seurat_obj$celltype <- Idents(seurat_obj)

pdf('umap_annotation.pdf')
DimPlot(seurat_obj, reduction = 'umap', label = TRUE, pt.size = 0.5) 
dev.off()

saveRDS(seurat_obj, file = 'seurat_aggr_ann.rds')

