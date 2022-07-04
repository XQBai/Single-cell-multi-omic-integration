library(Seurat)
library(factoextra)
library(dplyr)

seu <- readRDS('seurat_epi.rds')
subclone <- as.integer(seu$corrlabel)
RNAscale <- read.table('./cond/expr.infercnv.preliminary.dat', header = T)

# RNAscale <- seu@assays$RNA@data
chr_genes <- read.table('/mnt/ix1/Projects/M070_200622_GI_multiomics/scRNA/00_Code/GSVA_tables/Chromo_pathway.csv', header=T, sep=",", stringsAsFactors=F)

## Find the common genes with signal genes
signal_genes <- chr_genes %>% filter(Pathway %in% c('Chr3', 'Chr7', 'Chr8', 'Chr21'))
common_genes <- intersect(rownames(RNAscale), signal_genes$Genes)
RNA_mat <- RNAscale[common_genes, ]

d <- as.matrix(dist(t(RNA_mat)))
sil <- cluster::silhouette(subclone, dmatrix = d)
summary(sil)

p <- fviz_silhouette(sil)
pdf('silhoutte_plot.pdf')
p
dev.off()

cls.scatt <- clv::cls.scatt.data(t(RNAscale), subclone, dist = 'euclidean')
interclust = c("single", "complete", "average","centroid", "aveToCent", "hausdorff")
intraclust = c("complete","average","centroid")
davies1 <- clv.Davies.Bouldin(cls.scatt, intraclust, interclust)




