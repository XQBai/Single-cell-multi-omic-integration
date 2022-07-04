library(Seurat)
library(data.table)
#library(dbscan)
library("scales")
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
#setwd('/mnt/ix1/Projects/M070_200622_GI_multiomics/scDNA/A09_210112_CIN')
project_path <- '/mnt/ix1/Projects/M070_200622_GI_multiomics/scDNA'
seurat_path <- '/mnt/ix1/Projects/M070_200622_GI_multiomics/scDNA/A12_210825_seurat/'

# description of node bed file
# https://support.10xgenomics.com/single-cell-dna/software/pipelines/latest/output/bed
metrics_file <- paste0(project_path, '/A00_cellranger_dna/', opt$sample, '/outs/per_cell_summary_metrics.csv')
cellmetrics <- read.csv(metrics_file, sep=",", header=T)
nodebed <- read.table(paste0(project_path, '/A00_cellranger_dna/', opt$sample, '/outs/node_cnv_calls.bed'))
load(paste0(seurat_path, opt$sample, '/seurat_scDNA.Robj'))

colnames(nodebed) <- c('chrom', 'start', 'end', 'id', 'cnv', 'conf')
nodebed <- nodebed %>% filter(id <= max(id)/2)

singlenodebed <- nodebed %>% dplyr::filter(chrom != 'chrX' & chrom != 'chrY') %>% dplyr::filter(cnv!=2) %>% 
  dplyr::filter(id %in% which(seurat_scDNA$cellcomponment !='Noise')) %>% dplyr::group_by(id) %>% 
  dplyr::mutate(eventSize = (end-start)/1000) 

singlenodebed$amp.del <- matrix(0, nrow = dim(singlenodebed)[1], ncol = 1)
singlenodebed$amp.del[which(singlenodebed$cnv > 2)] <- 'amp'
singlenodebed$amp.del[which(singlenodebed$cnv < 2)] <- 'del'
singlenodebed$amp.del <- as.factor(singlenodebed$amp.del)

data <- singlenodebed
labels_order <- c('<=20kb', '<100kb', '<500kb', '<1Mb', '<5Mb', '<10Mb', '>=10Mb')

data[, c('chrom', 'start', 'end', 'id', 'cnv', 'conf')] <- NULL
ampdata <- data %>% filter(amp.del == 'amp')
ampFreq <- c()
ampFreq$freq <- hist(ampdata$eventSize, breaks = c(0, 20, 100, 500, 1000, 5000, 10000, max(ampdata$eventSize)))$counts
ampFreq$rate_cut <- labels_order
ampFreq$amp.del <- rep('amp', length(ampFreq$freq))
ampFreq <- as.data.frame(ampFreq)

deldata <- data %>% filter(amp.del == 'del')
delFreq <- c()
delFreq$freq <- hist(deldata$eventSize, breaks = c(0, 20, 100, 500, 1000, 5000, 10000, max(deldata$eventSize)))$counts
delFreq$rate_cut <- labels_order
delFreq$amp.del <- rep('del', length(delFreq$freq))
delFreq <- as.data.frame(delFreq)

dataFreq <- rbind(ampFreq, delFreq)
dataFreq <- dataFreq %>% filter(freq != 0)

dataFreq <- dataFreq %>% mutate(sumFreq=sum(freq)) %>% mutate(rate = round(freq/sumFreq,digits=4)*100)

#set order #########################################################################################
dataFreq$rate_cut <- factor(dataFreq$rate_cut, levels=labels_order,ordered = TRUE)
dataFreq <- dataFreq[order(dataFreq$freq, dataFreq$amp.del),]

#set plot text
plot_legend <- c("AMP", "DEL")
plot_title <- paste0("Distribution of Event Size")
text.offset <- max(round(dataFreq$rate,digits = 0))/25

mytheme <- theme_classic() + 
  theme(
    panel.background = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "gray95"),
    axis.ticks = element_blank(),
    text = element_text(family = "sans"),
    axis.title = element_text(color = "black", size = 12, face = 'bold'),
    axis.text = element_text(size = 10, color = "black", face = 'bold'),
    plot.title = element_text(size = 14, hjust = .5, color = "black", face = 'bold'),
    axis.line.y = element_line(size=1, linetype = 'dotted'),
    axis.line.x = element_blank(),
    axis.text.x = element_text(vjust = 0),
    axis.title.x = element_blank(),
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
    legend.position = c(0.1, 0.9),
    legend.text = element_text(color = "black", face = 'bold')
  )

p <- ggplot(dataFreq, aes(x=rate_cut, y=rate)) + 
  geom_bar(stat = "identity", aes(fill = amp.del), position = "dodge",width = 0.5) +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values = c("#00188F","#00BCF2")
                    ,breaks = c("amp","del")
                    ,labels = plot_legend
                    ,name = "") +
  geom_text(data = dataFreq
            , aes(label = comma(rate), y = rate +text.offset, color = amp.del)
            ,position = position_dodge(width =1)
            , size = 5) +
  scale_color_manual(values = c("#00BCF2", "#00188F"), guide = FALSE) +
  labs(y=" Proportation of event size (%)", title=plot_title) +
  #labs(x="Event Size",y="Percentage %", title=plot_title) +
  mytheme 
  # coord_flip() # change to Horizontal direction 
saveRDS(dataFreq, file = paste0(opt$sample, '_freq.rds'))
ggsave(paste0(opt$sample,'_event_size.pdf'), plot=p, dpi=300, units="in", height=7, width=7)


  