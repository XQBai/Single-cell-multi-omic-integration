library(pheatmap)

seu <- readRDS('seurat_epi.rds')
Idents(seu) <- seu@meta.data$corrlabel

hallmark_test <- read.csv("hallmark_test.csv")
pq_hallmark_df <- readRDS('pq_hallmark.rds')

hallmark_scaled <- GSVA_scaled(hallmark_test, seu)

#plotting
pdf("hallmark_gmt.pdf")
pheatmap::pheatmap(hallmark_scaled, fontsize =20, 
                   color=viridis(20),
                   fontsize_col = 15,
                   fontsize_row = 8
                   #cluster_rows = FALSE, cluster_cols = FALSE
)
dev.off()

#grabbing top_n to plot as above 
top3 <- rownames_to_column(hallmark_scaled, var = "pathway") %>% reshape2::melt()%>% group_by(variable)%>% top_n(3, value)
top3_hallmark_scaled <- as.data.frame(hallmark_scaled) %>% rownames_to_column(var="pathway")%>% dplyr::filter(pathway %in% top3$pathway)%>% column_to_rownames(var = "pathway")

pdf("hallmark_gmt_top3.pdf")
pheatmap::pheatmap(top3_hallmark_scaled, fontsize =15, 
                   color=viridis(20)
                   #cluster_rows = FALSE, cluster_cols = FALSE
)
dev.off()
