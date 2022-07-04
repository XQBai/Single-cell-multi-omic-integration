library(Seurat)

load('./P5931_tumor/seurat_scDNA.Robj')

AlleloscopeIden <- readRDS('./P5931_tumor/Alleloscope/rds/cell_lineage_identity.rds')

remain_cells <- intersect(which(seurat_scDNA$subclone != 'S phase') , 
                          intersect(which(seurat_scDNA$subclone != 'technical noise'), 
                                    which(seurat_scDNA$subclone != 'cellranger noise')))
seurat_scDNA_subset <- subset(seurat_scDNA, cells = remain_cells)

length(intersect(names(AlleloscopeIden), colnames(seurat_scDNA_subset)))

commoncells <- intersect(names(AlleloscopeIden), colnames(seurat_scDNA_subset))
seurat_common <- subset(seurat_scDNA_subset, cells = commoncells)
commonlabel <- as.data.frame(seurat_common$subclone)
commonlabel$allelo <- AlleloscopeIden[commoncells]
