library(dplyr)
library(pheatmap)

mean_dna_df <- readRDS('mean_dna_df.rds')
mean_rna_df <- readRDS('mean_rna_infercnv.rds')

mean_rna_df <- mean_rna_df[c(-1), ]
mean_dna_df <- mean_dna_df[rownames(mean_rna_df), ]

mean_dna_df <- round(mean_dna_df)

Similarity_mat <- matrix(0, nrow = dim(mean_dna_df)[1], ncol = dim(mean_rna_df)[1])
pvalue <- matrix(0, nrow = dim(mean_dna_df)[1], ncol = dim(mean_rna_df)[1])

for (i in 1:dim(Similarity_mat)[1]){
  for (j in 1:dim(Similarity_mat)[2]){
    Similarity_mat[i, j] <- cor(mean_dna_df[i, ], mean_rna_df[j, ], method = 'spearman')
    pvalue[i, j] <- cor.test(mean_dna_df[i, ], mean_rna_df[j, ], method = 'spearman')$p.value
  }
}

rownames(Similarity_mat) <- rownames(mean_dna_df)
colnames(Similarity_mat) <- rownames(mean_dna_df)

rownames(pvalue) <- rownames(mean_dna_df)
colnames(pvalue) <- rownames(mean_dna_df)

pdf('Similarity_metric.pdf')
pheatmap::pheatmap(Similarity_mat, cluster_rows = F, cluster_cols = F,
                   # color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
                   display_numbers = T, number_color = 'black',
                   main='Spearman Correlation', fontsize = 15,
                   fontsize_number = 15, fontsize_col = 15, fontsize_row = 15)
dev.off()

write.csv(pvalue, file = 'spearman_pvalue.csv')

