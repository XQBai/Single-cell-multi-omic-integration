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
library(RColorBrewer)
#修改legend https://www.shenxt.info/zh/post/ggplot2-square-legend/

setwd('/mnt/ix1/Projects/M070_200622_GI_multiomics/scDNA/A13_210910_mainfigure1/A01_210910_eventsize')

files <- Sys.glob("*_freq.rds")
samples <- gsub('_freq.rds', '', files)

# Merge frequency dataframe for all samples
alldata <- do.call("rbind", lapply(1:length(files), function(X){
  file <- readRDS(files[X])
  file$sample <-as.factor(rep(samples[X], dim(file)[1]))
  return(file)
  }))

levels(alldata$amp.del)[levels(alldata$amp.del) == 'amp'] <- 'Amp'
levels(alldata$amp.del)[levels(alldata$amp.del) == 'del'] <- 'Del'
names(alldata)[names(alldata) == 'sample'] <- 'Samples'

p <- ggplot(alldata, aes(rate_cut, rate, fill=Samples)) + 
  geom_bar(stat="identity", aes(fill = Samples)) +
  facet_grid(cols = vars(amp.del), scales = 'free_y') +
  ylab("Proportion of event size (%)") +
  labs(title = 'Comparison of CNV event size across samples') + 
  theme(plot.title = element_text(size = 20, color = 'black', face = 'bold', hjust = 0.5),
        axis.title.y = element_text(size = 20, color = 'black', face = 'bold', hjust = 0.5), 
        axis.text.y = element_text(size = 15, color = 'black', face = 'bold', hjust = 0.5),
        strip.text = element_text(size = 20, color = 'black', face = 'bold', hjust = 0.5),
        legend.title = element_text(size = 15, color = 'black', face = 'bold', hjust = 0.5), 
        legend.text = element_text(size = 15, color = 'black', face = 'bold', hjust = 0.5),
        axis.text.x = element_text(size = 15, angle = 45, color = 'black', face = 'bold', hjust = 1), 
        axis.title.x = element_blank() 
          )
p
ggsave(file= 'Samples_events_size.pdf', plot=p, dpi=300, units="in", height=8, width=16)

###################
# Plot the proportion of different event sizes of all samples
# Rename samples in alldata (dataframe)
alldata$Samples <- gsub('_tumor', '', alldata$Samples)
alldata$Samples[which(alldata$Samples == 'P5846')] <- 'P5846_gastric'
alldata$Samples[which(alldata$Samples == 'P5847')] <- 'P5847_gastric'
alldata$Samples[which(alldata$Samples == 'P5931')] <- 'P5931_gastric'
alldata$Samples[which(alldata$Samples == 'P6342')] <- 'P6342_gastric'

alldata$Samples[which(alldata$Samples == 'P5915_liver')] <- 'P5915_liver_met'
alldata$Samples[which(alldata$Samples == 'P6198_liver')] <- 'P6198_liver_met'
alldata$Samples[which(alldata$Samples == 'P6335_omentum')] <- 'P6335_omentum_met'

# Give an order of samples 
alldata$Samples <- factor(alldata$Samples, levels = c('P5846_gastric', 'P5847_gastric', 'P5931_gastric',
                                                      'P6342_gastric', 'P6461_colon', 'P6461_liver', 'P6593_colon', 'P6593_liver_met',
                                                      'P5915_liver_met', 'P6198_liver_met', 'P6335_omentum_met'))

levels(alldata$rate_cut)[levels(alldata$rate_cut) == '<=20kb'] <- '[0, 20kb)'
levels(alldata$rate_cut)[levels(alldata$rate_cut) == '<100kb'] <- '[20kb, 100kb)'
levels(alldata$rate_cut)[levels(alldata$rate_cut) == '<500kb'] <- '[100kb, 500kb)'
levels(alldata$rate_cut)[levels(alldata$rate_cut) == '<1Mb'] <- '[500kb, 1Mb)'
levels(alldata$rate_cut)[levels(alldata$rate_cut) == '<5Mb'] <- '[1Mb, 5Mb)'
levels(alldata$rate_cut)[levels(alldata$rate_cut) == '<10Mb'] <- '[5Mb, 10Mb)'
levels(alldata$rate_cut)[levels(alldata$rate_cut) == '>=10Mb'] <- '[10Mb, Inf)'



alldata_int <- alldata %>% dplyr::group_by(Samples, rate_cut) %>% dplyr::mutate(rate_int = sum(rate)) %>% 
  dplyr::filter(!duplicated(rate_int))

p1 <- ggplot() + 
  geom_bar(aes(x = Samples, y = rate_int, fill = rate_cut), data = alldata_int,
           colour = 'black', stat = 'identity', position = 'stack') + 
  scale_fill_manual(values = brewer.pal(7, 'Set2'))  + 
  ylab("Proportion (%)") +
  labs(title = 'Various event sizes with aneuploidy across samples') + 
  theme(plot.title = element_text(size = 18, color = 'black', face = 'bold', hjust = 0.5),
        axis.title.y = element_text(size = 20, color = 'black', face = 'bold', hjust = 0.5), 
        axis.text.y = element_text(size = 15, color = 'black', face = 'bold', hjust = 0.5),
        strip.text = element_text(size = 20, color = 'black', face = 'bold', hjust = 0.5),
        legend.title = element_text(size = 15, color = 'black', face = 'bold', hjust = 0.5), 
        legend.text = element_text(size = 15, color = 'black', face = 'bold', hjust = 0.5),
        axis.text.x = element_text(size = 15, angle = 45, color = 'black', face = 'bold', hjust = 1), 
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        legend.key.size  = unit(2, 'line'),
      
  ) + 
  guides(fill = guide_legend(title = "Event sizes")) 

p1
ggsave(file= 'All_Samples_events_size.pdf', plot=p1, dpi=300, units="in", height=7, width=10)

##### Plot the proportion of amplification and deletions 
alldata_ampdel <- alldata %>% dplyr::group_by(Samples, amp.del) %>% dplyr::mutate(rate_int = sum(rate)) %>%
  dplyr::filter(!duplicated(amp.del))
alldata_ampdel$rate_int <- round(alldata_ampdel$rate_int, 1)

text.offset <- max(round(alldata_ampdel$rate_int,digits = 0))/50

p2 <- ggplot(alldata_ampdel, aes(x=Samples, y=rate_int)) + 
  geom_bar(stat = "identity", aes(fill = amp.del), position = "dodge") +
  scale_fill_manual(values = brewer.pal(7, 'Set1'), breaks = c("Amp","Del"))  +
  ylab("Proportion (%)") +
  labs(title = 'Comparison of aneuploidy events across samples') + 
  geom_text(data = alldata_ampdel
            , aes(label = comma(rate_int, accuracy = 0.1), y = rate_int + text.offset, color = amp.del)
            , show.legend = FALSE, position = position_dodge(width = 1)
            , size = 5) + 
  theme(plot.title = element_text(size = 18, color = 'black', face = 'bold', hjust = 0.5),
        axis.title.y = element_text(size = 20, color = 'black', face = 'bold', hjust = 0.5),
        axis.text.y = element_text(size = 15, color = 'black', face = 'bold', hjust = 0.5),
        axis.line.y = element_line(linetype = 1, color = 'black', size = 0.5), 
        strip.text = element_text(size = 20, color = 'black', face = 'bold', hjust = 0.5),
        legend.title = element_text(size = 20, color = 'black', face = 'bold', hjust = 0.5),
        legend.text = element_text(size = 15, color = 'black', face = 'bold', hjust = 0.5),
        axis.text.x = element_text(size = 15, angle = 45, color = 'black', face = 'bold', hjust = 1),
        axis.title.x = element_blank(),
        axis.line.x = element_line(linetype = 1, color = 'black', size = 0.5), 
        panel.background = element_blank(),
        legend.key.size  = unit(2, 'line')
        ) + 
  guides(fill = guide_legend(title = "Status"))
p2
ggsave(file= 'All_samples_amp_del.pdf', plot=p2, dpi=300, units="in", height=6, width=8)

## Combine two figures
p1 + p2
ggsave(file= 'combined_event_size_samples.pdf', plot=p1 + p2, dpi=300, units="in", height=8, width=20)


