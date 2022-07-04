library(ape)
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(ggplot2)


load('seurat_scDNA.Robj')
load('seurat_scDNA_tumor.Robj')

labels <- as.character(Idents(seurat_scDNA_subset))
seurat_scDNA_subset@meta.data$subclones <- labels

# label <- as.numeric(Idents(seurat_scDNA_subset)) 
# seurat_scDNA_subset@meta.data$subclones <- paste0('C', label)
# labels <- seurat_scDNA_subset$subclones

if(length(unique(labels)) > 9){
  mycol = brewer.pal(9, 'Set1')
  Label_color <- c(mycol, brewer.pal(12, 'Set3')[10:12], "#666666", brewer.pal(5, 'Set1'))
  names(Label_color) <- sort(unique(labels))
}else if(length(unique(labels)) >= 3 & length(unique(labels)) <= 9){
  Label_color <- brewer.pal(length(unique(labels)), "Set1")
  names(Label_color) <- sort(unique(labels))
}else if(length(unique(labels)) == 2){
  Label_color <- c("#e41a1c", "#377eb8")
  names(Label_color) <- sort(unique(labels))
}else if(length(unique(labels)) == 1){
  Label_color <- c('#e41a1c')
  names(Label_color) <- sort(unique(labels))
}

pdf('scDNA_umap_no_normal.pdf')
DimPlot(seurat_scDNA_subset, group.by = 'subclones', cols = Label_color, reduction = "umap") +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  labs(title = '') + 
  guides(color=guide_legend("Subclone")) +
  theme(plot.title = element_text(size = 11, color = 'black', face = 'bold', hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 15, color = 'black', face = 'bold')) +
  theme(axis.text.y = element_text(size = 15, color = 'black', face = 'bold')) +
  theme(axis.title.x = element_text(size = 15, color = 'black', face = 'bold')) + 
  theme(axis.title.y = element_text(size = 15, color = 'black', face = 'bold')) + 
  theme(legend.text= element_text(size=15,color="black", face= "bold", vjust=0.5, hjust=0.5)) +
  theme(legend.title = element_text(size=15,color="black", face= "bold"))
dev.off()


# # labels <- seurat_scDNA$subclone
# label <- as.numeric(Idents(seurat_scDNA))
# seurat_scDNA@meta.data$subclones <- paste0('C', label)
# labels <- seurat_scDNA$subclones
# 
# if(length(unique(labels)) > 9){
#   mycol = brewer.pal(9, 'Set1')
#   Label_color <- c(mycol, brewer.pal(12, 'Set3')[10:12], "#666666")
#   names(Label_color) <- sort(unique(labels))
# }else if(length(unique(labels)) >= 3 & length(unique(labels)) <= 9){
#   Label_color <- brewer.pal(length(unique(labels)), "Set1")
#   names(Label_color) <- sort(unique(labels))
# }else if(length(unique(labels)) == 2){
#   Label_color <- c("#e41a1c", "#377eb8")
#   names(Label_color) <- sort(unique(labels))
# }else if(length(unique(labels)) == 1){
#   Label_color <- c('#e41a1c')
#   names(Label_color) <- sort(unique(labels))
# }
# 
# pdf('scDNA_umap_all.pdf')
# DimPlot(seurat_scDNA, group.by = 'subclone', cols = Label_color, reduction = "umap") +
#   xlab("UMAP 1") + 
#   ylab("UMAP 2") +
#   labs(title = '') + 
#   guides(color=guide_legend("Subclones")) +
#   theme(plot.title = element_text(size = 11, color = 'black', face = 'bold', hjust = 0.5)) +
#   theme(axis.text.x = element_text(size = 15, color = 'black', face = 'bold')) +
#   theme(axis.text.y = element_text(size = 15, color = 'black', face = 'bold')) +
#   theme(axis.title.x = element_text(size = 15, color = 'black', face = 'bold')) + 
#   theme(axis.title.y = element_text(size = 15, color = 'black', face = 'bold')) + 
#   theme(legend.text= element_text(size=15,color="black", face= "bold", vjust=0.5, hjust=0.5)) +
#   theme(legend.title = element_text(size=15,color="black", face= "bold"))
# dev.off()

