# run infercnv when use one of clutser in seurat object as normal cells (reference)
library(Seurat)
library(optparse)
library(dplyr)
library(viridis)
library(cowplot)
library(future)
library(infercnv)

options(future.globals.maxSize = 6000 * 1024^2)
plan('multiprocess', workers= 18)

seurat_obj <- readRDS('seurat_epi.rds')

# a raw counts matrix of single-cell RNA-Seq expression
DefaultAssay(seurat_obj) <- 'RNA'
counts_matrix = GetAssayData(seurat_obj, slot='counts')

# create the infercnv object 
infercnv_obj = CreateInfercnvObject(counts_matrix, 
                                    annotations_file= "idents.txt", 
                                    delim="\t", gene_order_file="gene_locs.sorted.bed", 
                                    ref_group_names=NULL)

# perform infercnv operations to reveal cnv signal
infercnv_obj_final = infercnv::run(infercnv_obj, 
                                   cutoff=.1, 
                                   out_dir="test", 
                                   cluster_by_groups=T, 
                                   #num_threads=40, 
                                   denoise=T,
                                   HMM=T)
