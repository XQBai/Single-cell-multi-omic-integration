library(dplyr)
library(Seurat)
library(ggplot2)

seurat <- readRDS('seurat_epi.rds')
dna <- readRDS('mean_dna_df.rds')
dna <- round(dna)
label <- as.character(seurat$corrlabel)
label[which(label == 'normal')] <- 'Normal'
rna_df <- readRDS('rna_infercnv_mean_chrom.rds')
rna <- t(rna_df)

rsme <- function(rna, label, dna){
  
  rmsd <- matrix(0, nrow = length(unique(label)))
  for (i in 1:length(unique(label))){
    
    s <- as.character(unique(label)[i])
    index <- which(label == s)
    sc_rna <- rna[, index]
    sub_dna <- dna[s, ]
    
    ones <- matrix(1, nrow = 1, ncol = dim(sc_rna)[2])
    diff_mat <- as.matrix(sc_rna - as.matrix(sub_dna) %*% ones)
    l2_norm <- apply(diff_mat, 2, function(x){norm(x, type = '2')})
    rmsd[i] <- sqrt(sum(l2_norm)/length(index))
    # rmsd[i] <- sqrt(sum(l2_norm))
  }
  rmsd_avg <- mean(rmsd)
  
  return(rmsd_avg)
}

sig_chr <- c('chr3', 'chr7', 'chr8', 'chr20', 'chr21')
rna_signal <- rna[sig_chr, ]
dna_signal <- dna[, sig_chr]

permutation_rmsd <- matrix(0, nrow = 1000)

for (i in 1:1000){
  Label <- sample(unique(label), length(label), replace = T)
  rmsd_avg <- rsme(rna_signal, Label, dna_signal)
  permutation_rmsd[i] <- rmsd_avg
}

write.csv(permutation_rmsd, file = 'permutation_rmsd.csv')

pdf('permutation_density.pdf')
plot(density(permutation_rmsd))
dev.off()

### Plot only for epithelial cells 
per_rsmd <- read.csv('permutation_rmsd.csv', header = F, sep = ',', 
                     skip = 1)
rsmd_1000 <- per_rsmd$V2
rsmd_obv <- 1.549
z <- (rsmd_obv - mean(rsmd_1000))/(sd(rsmd_1000)/sqrt(1000))
pvalue <- 2 * pnorm(-abs(z))

df <- data.frame(x = rsmd_1000)

p <- ggplot(df, aes(x = x)) + 
  geom_histogram(color = 'black', fill = '#a6bddb', bins = 50) + 
  geom_vline(aes(xintercept = rsmd_obv), color = '#fc9272', 
             linetype = 'dashed', size = 0.5) + 
  theme_classic() + 
  # xlim(c(1.54, 1.61)) +
  ylab('count') + 
  xlab('Root mean square error') + 
  theme(legend.position = 'top')
p

pdf('permutation.pdf')
p
dev.off()
