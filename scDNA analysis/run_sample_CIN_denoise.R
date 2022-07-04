library(Seurat)
library(data.table)
#library(dbscan)
library(dplyr)
library(Matrix)
library(mixtools)
library(optparse)
library(ggplot2)
library(future)

options(future.globals.maxSize = 24000 * 1024^2)
plan('multiprocess', workers=24)

option_list <- list(
  make_option(c('-s', '--sample'), action='store', default = NULL, type='character', help = 'Sample name')
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$sample)) {
  print_help(opt_parser)
  stop("Need sample name.")
}

print(opt)
project_path <- '/mnt/ix1/Projects/M070_200622_GI_multiomics/scDNA'
seurat_path <- '/mnt/ix1/Projects/M070_200622_GI_multiomics/scDNA/A12_210825_seurat/'

metrics_file <- paste0(project_path, '/A00_cellranger_dna/', opt$sample, '/outs/per_cell_summary_metrics.csv')
cellmetrics <- read.csv(metrics_file, sep=",", header=T)
nodebed <- read.table(paste0(project_path, '/A00_cellranger_dna/', opt$sample, '/outs/node_cnv_calls.bed'))
load(paste0(seurat_path, opt$sample, '/seurat_scDNA.Robj'))

colnames(nodebed) <- c('chrom', 'start', 'end', 'id', 'cnv', 'conf')
nodebed <- nodebed %>% filter(id <= max(id)/2)

singlenodebed <- nodebed %>% dplyr::group_by(id) %>% dplyr::mutate(segconut = length(id)) %>% 
  dplyr::mutate(CINcount=length(which(cnv!=2)==TRUE)) %>%
  dplyr::mutate(CIN_percent = CINcount/segconut) %>%
  dplyr::filter(!duplicated(id)) %>% arrange(id) 

singlenodebed$barcode <- cellmetrics$barcode
singlenodebed$Noise <- seurat_scDNA$cellcomponment
# singlenodebed$Noise[which(singlenodebed$Noise != 'Noise')] <- 'No'
# singlenodebed$Noise[which(singlenodebed$Noise == 'Noise')] <- 'Yes'
# 
# p <- ggplot(data = singlenodebed, mapping = aes(x = CIN_percent)) + 
#   geom_density(colour="steelblue", fill='lightblue') + 
#   xlab('CIN Percent') + 
#   ylab("Density") + 
#   theme(plot.title = element_text(size = 11, color = 'black', face = 'bold', hjust = 0.5)) +
#   theme(axis.text.x = element_text(size = 11, color = 'black', face = 'bold')) +
#   theme(axis.text.y = element_text(size = 11, color = 'black', face = 'bold')) +
#   theme(axis.title.x = element_text(size = 15, color = 'black', face = 'bold')) + 
#   theme(axis.title.y = element_text(size = 15, color = 'black', face = 'bold')) + 
#   theme(legend.text= element_text(size=15,color="black", face= "bold", vjust=0.5, hjust=0.5)) +
#   theme(legend.title = element_text(size=15,color="black", face= "bold"))
# 
# p1 <- ggplot(data = singlenodebed, mapping = aes(x = CIN_percent, fill=Noise)) + 
#   geom_density(alpha=0.5) + 
#   xlab('CIN Percent') + 
#   ylab("Density") + 
#   theme(plot.title = element_text(size = 11, color = 'black', face = 'bold', hjust = 0.5)) +
#   theme(axis.text.x = element_text(size = 11, color = 'black', face = 'bold')) +
#   theme(axis.text.y = element_text(size = 11, color = 'black', face = 'bold')) +
#   theme(axis.title.x = element_text(size = 15, color = 'black', face = 'bold')) + 
#   theme(axis.title.y = element_text(size = 15, color = 'black', face = 'bold')) + 
#   theme(legend.text= element_text(size=15,color="black", face= "bold", vjust=0.5, hjust=0.5)) +
#   theme(legend.title = element_text(size=15,color="black", face= "bold"))

# saveRDS(nodebed, file = paste0(opt$sample, '_nodebed.rds'))
saveRDS(singlenodebed, file = paste0(opt$sample,'_singlenodebed.rds' ))
# ggsave(paste0(opt$sample, '_CIN_fraction.png'), plot=p, dpi=300, units="in", height=6, width=8)
# ggsave(paste0(opt$sample, '_CIN_fraction_noise.png'), plot=p1, dpi=300, units="in", height=6, width=8)


