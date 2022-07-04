library(pheatmap)
library(ggplot2)
library(Seurat)
library(stringr)
library(Rtsne)
#library(uwot)
library(mclust)
library(optparse)
library(biomaRt)
library(cowplot)
library(ggplotify)
library(RColorBrewer)
library(minerva)
library(data.table)
library(biomaRt)
library(future)
# library(GRanges)
library(GenomicRanges)
source('/home/xiangqi/shared/R_scripts/InputData.R')
source('/home/xiangqi/shared/R_scripts/InitializeCNV.R')

options(future.globals.maxSize = 10000*1024*2)
plan('multiprocess', workers = 24)

option_list <- list(
  make_option(c("-i", "--input"), action="store", default=NULL, type='character', help="Input base folder for mtx and tsv files"),
  make_option(c("-s", "--sample"), action="store", default=NULL, type='character', help="Sample name/output prefix")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Need input folder.")
}

if (is.null(opt$sample)) {
  print_help(opt_parser)
  stop("Need sample name.")
}

cnvpos <- read.table('GRCh38_cellranger_20k.canonical.rownumbers.bed')
names(cnvpos) <- c('chr', 'start', 'end', 'bin')

armpos <- read.table('/mnt/ix1/Projects/M070_200622_GI_multiomics/Integrated/TCGA_int/chromosome_arm_positions_grch38.txt', header = TRUE)
names(armpos) <- c('idf', 'chr', 'start', 'end', 'length')
armpos$chr <- paste0('chr', armpos$chr)

DNAmatrix <- fread(opt$input)

## Remove the cnv on sex chromosome
cnvpos <- cnvpos %>% filter(chr != 'chrX' & chr != 'chrY')
cnvpos$bin <- cnvpos$bin + 1

## Convert these to GRanges object
cnvgr <- makeGRangesFromDataFrame(cnvpos, keep.extra.columns = TRUE)
armgr <- makeGRangesFromDataFrame(armpos, keep.extra.columns = TRUE)

## Find the overlaps between cnv position and chr arm position
olaps <- findOverlaps(cnvgr, armgr)

# the corresponding arm for cnv
armbin <- data_frame(bin = cnvpos$bin[queryHits(olaps)],
                      arm = armpos$idf[subjectHits(olaps)])

DNAmatrix <- DNAmatrix[armbin$bin, ]
DNAmatrix$arm <- armbin$arm

DNAmatrixArm <- DNAmatrix %>% group_by(arm) %>% summarise_all(list(mean), na.rm = T)
## Rearrange the arm according to order
chr <- gsub('p', '', DNAmatrixArm$arm)
chr <- gsub('q', '', chr)
DNAmatrixArm$chr <- as.numeric(chr)
DNAmatrixArm <- DNAmatrixArm %>% arrange(chr)
DNAmatrixArm$arm <- NULL
DNAmatrixArm$chr <- NULL
DNAmatrixArm <- DNAmatrixArm %>% summarise_all(round)

save(DNAmatrixArm, file=paste0(opt$sample,'_CNV_Arm_matrix.Robj'))





