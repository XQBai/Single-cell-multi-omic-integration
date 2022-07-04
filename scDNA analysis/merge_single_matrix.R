library(dplyr)

setwd('/mnt/ix1/Projects/M070_200622_GI_multiomics/scDNA/A13_210910_mainfigure1/A04_210910_TCGA/A01_210910_armMatrix')

#### Primary Single cells samples
load('./P5846_tumor_CNV_Arm_matrix_nonoise.Robj')
P5846_matrix <- DNAmatrixArm

load('./P5847_tumor_CNV_Arm_matrix_nonoise.Robj')
P5847_matrix <- DNAmatrixArm

load('./P5931_tumor_CNV_Arm_matrix_nonoise.Robj')
P5931_matrix <- DNAmatrixArm

load('./P6342_tumor_CNV_Arm_matrix_nonoise.Robj')
P6342_matrix <- DNAmatrixArm

# load('./P6461_tumor_colon_CNV_Arm_matrix_nonoise.Robj')
# P6461_tumor_colon_matrix <- DNAmatrixArm
# 
# load('./P6593_colon_tumor_CNV_Arm_matrix_nonoise.Robj')
# P6593_colon_tumor_matrix <- DNAmatrixArm

SinglecellArm_matrix <- cbind(cbind(cbind(cbind(cbind(P5846_matrix, P5847_matrix), P5931_matrix), P6342_matrix)))
armpos <- read.table('/mnt/ix1/Projects/M070_200622_GI_multiomics/Integrated/TCGA_int/chromosome_arm_positions_grch38.txt', header = TRUE)
rownames(SinglecellArm_matrix) <- armpos$Idf
colname <- rep(c('P5846', 'P5847', 'P5931', 'P6342'), c(dim(P5846_matrix)[2], dim(P5847_matrix)[2], 
                                                                                      dim(P5931_matrix)[2], dim(P6342_matrix)[2]))

colnames(SinglecellArm_matrix) <- paste0(colname, '_cell', 1:dim(SinglecellArm_matrix)[2])

saveRDS(SinglecellArm_matrix, file = './SinglecellArm_matrix.rds')


# ############# Gene matrix
# setwd('/mnt/ix1/Projects/M070_200622_GI_multiomics/Integrated/TCGA_int/A02_20210707_SingleCell_processed/A10_210707_geneMatrix')
# 
# #### Primary Single cells samples
# load('./P5846_tumor_CNV_gene_matrix_nonoise.Robj')
# P5846_matrix <- CNV_gene_Matrix
# 
# load('./P5847_tumor_CNV_gene_matrix_nonoise.Robj')
# P5847_matrix <- CNV_gene_Matrix
# 
# load('./P5931_tumor_CNV_gene_matrix_nonoise.Robj')
# P5931_matrix <- CNV_gene_Matrix
# 
# load('./P6342_tumor_CNV_gene_matrix_nonoise.Robj')
# P6342_matrix <- CNV_gene_Matrix
# 
# load('./P6461_tumor_colon_CNV_gene_matrix_nonoise.Robj')
# P6461_tumor_colon_matrix <- CNV_gene_Matrix
# 
# load('./P6593_colon_tumor_CNV_gene_matrix_nonoise.Robj')
# P6593_colon_tumor_matrix <- CNV_gene_Matrix
# 
# SinglecellGene_matrix <- cbind(cbind(cbind(cbind(cbind(P5846_matrix, P5847_matrix), P5931_matrix), P6342_matrix), P6461_tumor_colon_matrix), P6593_colon_tumor_matrix)
# 
# colname <- rep(c('P5846', 'P5847', 'P5931', 'P6342', 'P6461_colon', 'P6593_colon'), c(dim(P5846_matrix)[2], dim(P5847_matrix)[2], 
#                                                                                       dim(P5931_matrix)[2], dim(P6342_matrix)[2],
#                                                                                       dim(P6461_tumor_colon_matrix)[2], dim(P6593_colon_tumor_matrix)[2]))
# 
# colnames(SinglecellGene_matrix) <- paste0(colname, '_cell', 1:dim(SinglecellGene_matrix)[2])
# 
# genename <- read.table('genelist.txt')
# rownames(SinglecellGene_matrix) <- genename$V1
# saveRDS(SinglecellGene_matrix, file = './SinglecellGene_matrix.rds')




















