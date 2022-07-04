library(Seurat)
library(dplyr) 
library(Matrix)
library(GSEABase)
library(GSVA)

#load seurat object
seu <- readRDS("../HCT116_NGG_1.seurat_dfx.annot.subset.group.rds")
hallmark <- getGmt("/mnt/ix1/Resources/scRNA_Ref/GSEA_Genesets/h.all.v6.2.symbols.gmt")
HCT_NGG_1_hallmark <- gsva(as.matrix(GetAssayData(seu, slot = "data", assay = "RNA")), hallmark, kcdf="Gaussian",  mx.diff=T, verbose=FALSE, parallel.sz=3, min.sz=10)
saveRDS(HCT_NGG_1_hallmark, file="HCT_NGG_1_hallmark.rds", version =2)

