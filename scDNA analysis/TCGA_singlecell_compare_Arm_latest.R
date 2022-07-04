library(ggplot2)
library(reshape2)
library(ggthemes)
library(dplyr)
library(stringr)
library(cowplot)
library(GenomicRanges)
library(future)
library(pheatmap)

setwd('/mnt/ix1/Projects/M070_200622_GI_multiomics/scDNA/A13_210910_mainfigure1/A04_210910_TCGA')

tcgaArm <- read.csv('./stad_tcga_combine_arm.csv')
chrArm <- tcgaArm$Cytoband
tcgaArm$X <- NULL
# tcgaArm$Cytoband <- NULL
rownames(tcgaArm) <- chrArm

# convert the wide matrix to long dataframe
tcga_df <- tcgaArm %>% tidyr::pivot_longer(cols = starts_with('TCGA'),
                                     names_to = 'Sample', values_to = 'cn')
# tcga_df$Cytoband <- as.factor(tcga_df$Cytoband)
tcga_df$Sample <- as.factor(tcga_df$Sample)

### Have to use .drop=FALSE if we want to use group information
tcga_ct <- tcga_df %>% dplyr::group_by(Cytoband, .drop = FALSE) %>% dplyr::mutate(Amp = sum(cn > 0), Del = sum(cn < 0), Normal = sum(cn == 0)) %>% ungroup()
tcga_ct <- tcga_ct %>% dplyr::filter(!duplicated(Cytoband))
tcga_ct <- tcga_ct[c('Cytoband', 'Amp', 'Del', "Normal")]

# Convert the wide data frame to long data frame
tcga_ct <- tcga_ct %>% tidyr::gather(Event, cn, Amp:Normal) 
tcga_ct <- tcga_ct %>% filter(Event != 'Normal')
# tcga_ct$order <- c(seq(1, 40, 1), seq(1, 40, 1))
tcga_ct$cn <- tcga_ct$cn/411 * 100
# tcga_ct$Event <- as.factor(tcga_ct$Event)
tcga_ct_top <- tcga_ct %>% dplyr::group_by(Event, .drop = FALSE) %>% top_n(10, cn) %>% ungroup()

# https://stackoverflow.com/questions/49290010/how-to-part-diverging-bar-plots-in-r
# you'll want to left-justify some labels and right-justify others, so make that a variable now:
tcga_ct_top$cn <- ifelse(tcga_ct_top$Event == 'Del', tcga_ct_top$cn*-1, tcga_ct_top$cn)
tcga_ct_top$just <- ifelse(tcga_ct_top$Event == 'Del', 0, 1)
tcga_ct_top <- tcga_ct_top[order(tcga_ct_top$cn), ]
tcga_ct_top$Event <- factor(tcga_ct_top$Event, levels = tcga_ct_top$Event)
# tcga_ct_top$Cytoband <- factor(tcga_ct_top$Cytoband, levels = tcga_ct_top$Cytoband)

p1 <- ggplot(tcga_ct_top, aes(x = Cytoband, y= cn, fill = Event)) + 
  geom_bar(stat='identity', width = .7) + 
  theme_bw() +
  labs(title = 'Chromosomal arm aneuploidy events of stomach adenocarcinoma samples in TCGA', x="",
       y="Percentage of samples (%)") + 
  scale_y_continuous(labels = c(40, 20, 0, 20, 40, 60), breaks=c(-40, -20, 0, 20, 40, 60)) +
  scale_fill_manual(values =  c('#fb8072', '#80b1d3')) + 
  theme(plot.title = element_text(size = 15, color = 'black', face = 'bold', hjust = 0.5),
        axis.title.x = element_text(size = 12, color = 'black', face = 'bold', hjust = 0.5),
        axis.text.x = element_text(size = 12, color = 'black', face = 'bold', hjust = 0.5),
        legend.text= element_text(size=12, color="black", face= "bold", vjust=0.5, hjust=0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank()
        ) +
  geom_text(aes(x = Cytoband, y = 0, label = Cytoband), size = 5, hjust=tcga_ct_top$just) + 
  coord_flip()

p1 


#### Single cells samples
scArm <- readRDS('./A01_210910_armMatrix/SinglecellArm_matrix.rds')
scArm$Cytoband <- rownames(scArm)
sc_df <- scArm %>% tidyr::pivot_longer(cols = starts_with('P'), names_to = 'Sample', values_to = 'cn')
# sc_df$Cytoband <- as.factor(sc_df$Cytoband)
sc_df$Sample <- as.factor(sc_df$Sample)
## Need to delete the 13p, 14p, 15p and 22p because they are too short
sc_df <- sc_df %>% dplyr::filter(Cytoband %in% tcgaArm$Cytoband)

### Have to use .drop=FALSE if we want to use group information
sc_ct <- sc_df %>% dplyr::group_by(Cytoband, .drop = FALSE) %>% dplyr::mutate(Amp = sum(cn > 2), Del = sum(cn < 2), Normal = sum(cn == 2)) %>% ungroup()
sc_ct <- sc_ct %>% dplyr::filter(!duplicated(Cytoband))
sc_ct <- sc_ct[c('Cytoband', 'Amp', 'Del', "Normal")]

# Convert the wide data frame to long data frame
sc_ct <- sc_ct %>% tidyr::gather(Event, cn, Amp:Normal) 
sc_ct <- sc_ct %>% dplyr::filter(Event != 'Normal')
# tcga_ct$order <- c(seq(1, 40, 1), seq(1, 40, 1))
sc_ct$cn <- sc_ct$cn/length(unique(sc_df$Sample)) * 100
# sc_ct$Event <- as.factor(sc_ct$Event)
sc_ct_top <- sc_ct %>% dplyr::group_by(Event, .drop = FALSE) %>% top_n(10, cn) %>% ungroup()

# https://stackoverflow.com/questions/49290010/how-to-part-diverging-bar-plots-in-r
# you'll want to left-justify some labels and right-justify others, so make that a variable now:
sc_ct_top$cn <- ifelse(sc_ct_top$Event == 'Del', sc_ct_top$cn*-1, sc_ct_top$cn)
sc_ct_top$just <- ifelse(sc_ct_top$Event == 'Del', 0, 1)
sc_ct_top <- sc_ct_top[order(sc_ct_top$cn), ]
# sc_ct_top$Event <- factor(sc_ct_top$Event, levels = sc_ct_top$Event)
sc_ct_top$Cytoband <- factor(sc_ct_top$Cytoband, levels = sc_ct_top$Cytoband)

p2 <- ggplot(sc_ct_top, aes(x = Cytoband, y= cn, fill = Event)) + 
  geom_bar(stat='identity', width = .7) + 
  theme_bw() +
  labs(title = 'Chromosomal arm aneuploidy events of primary gastric cancer single cells', x="",
       y="Percentage of cells (%)") + 
  scale_y_continuous(labels = c('100', '80', '60', '40', '20', '0', '20', '30'), breaks = c(-100, -80, -60, -40, -20, 0, 20, 30)) +
  scale_fill_manual(values =  c('#fb8072', '#80b1d3')) + 
  theme(plot.title = element_text(size = 15, color = 'black', face = 'bold', hjust = 0.5),
        axis.title.x = element_text(size = 12, color = 'black', face = 'bold', hjust = 0.5),
        axis.text.x = element_text(size = 12, color = 'black', face = 'bold', hjust = 0.5),
        legend.text= element_text(size=12, color="black", face= "bold", vjust=0.5, hjust=0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank()
  ) +
  geom_text(aes(x = Cytoband, y = 0, label = Cytoband), size = 5, hjust=sc_ct_top$just) +
  coord_flip()

p2
ggsave(filename = 'scdna.pdf', plot = p2, width = 8.5, height = 6, dpi = 300)
ggsave(filename = 'tcga.pdf', plot = p1, width = 9.3, height = 6, dpi = 300)


