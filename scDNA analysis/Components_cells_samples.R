library(ggplot2)
library(gdata)
library(dplyr)
library(readxl)
library(RColorBrewer)

setwd('/mnt/ix1/Projects/M070_200622_GI_multiomics/scDNA/A13_210910_mainfigure1/A03_210910_cellcomponment')
data <- read_excel('/mnt/ix1/Projects/M070_200622_GI_multiomics/scDNA/A13_210910_mainfigure1/A03_210910_cellcomponment/Supplementay_tables/Supplement_tables.xlsx',
                   sheet = 'Supplementary Table S4', skip = 1)

data$Sample <- c('P5846_gastric', 'P5847_gastric', 'P5931_gastric', 'P6342_gastric', 'P6461_colon', 'P6461_liver', 'P6593_colon', 'P6593_liver_met', 'P5915_liver_met',
                 'P6198_liver_met', 'P6335_omentum_met')
# data$freq_noise <- data$`#Noisy cells`/data$`#Total cells`
# data$freq_normal <- data$`#Normal cells`/data$`#Total cells`
# data$freq_tumor <- data$`#Pure tumor cells`/data$`#Total cells`

Sample <- c()
Fraction <- c()
for (i in 6:dim(data)[2]){
  Fraction <- c(Fraction, unlist(data[ ,i]))
  Sample <- c(Sample, as.character(data$Sample))
}
Type <- rep(c('Normal', 'G0/G1', 'S Phase'), c(length(data$Sample), length(data$Sample), length(data$Sample)))

alldata <- data.frame(Sample=Sample, Fraction = round(Fraction,2), Type=Type)
alldata$Sample <- factor(alldata$Sample, levels = c('P6335_omentum_met', 'P6198_liver_met', 'P5915_liver_met','P6593_liver_met', 
                                                    'P6593_colon', 'P6461_liver', 'P6461_colon', 'P6342_gastric', 'P5931_gastric', 
                                                    'P5847_gastric', 'P5846_gastric'))
alldata$Type <- factor(alldata$Type, levels = c('S Phase', 'G0/G1', 'Normal'))
alldata$Fraction <- alldata$Fraction * 100


p <- ggplot(alldata, mapping = aes(x = Sample, y = Fraction, fill = Type)) +
  geom_bar(stat = "identity",  width=0.5) +
  coord_flip() + 
  scale_fill_manual(values = c("#FC8D62", "#66C2A5", "#8DA0CB"), 
                    name = 'Cell Type') + 
  labs(y='Fraction (%)', title = 'Variable fractions of cell types') + 
  theme(plot.title = element_text(size = 20, color = 'black', face = 'bold', hjust = 0.5),
        axis.title.x = element_text(size = 18, color = 'black', face = 'bold', hjust = 0.5), 
        axis.text.y = element_text(size = 15, color = 'black', face = 'bold', hjust = 0.5),
        strip.text = element_text(size = 18, color = 'black', face = 'bold', hjust = 0.5),
        legend.title = element_text(size = 18, color = 'black', face = 'bold', hjust = 0.5), 
        legend.text = element_text(size = 18, color = 'black', face = 'bold', hjust = 0.5),
        legend.position = 'top', 
        axis.text.x = element_text(size = 15, angle = 45, color = 'black', face = 'bold', hjust = 1), 
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),#只删掉y轴刻度线
        axis.line.x = element_line(size = 0.5, colour = 'black'), #加上x轴的坐标轴
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) + 
  guides(fill = guide_legend(reverse = TRUE)) #颠倒图例的顺序

ggsave(file = 'Cell_fraction.pdf', plot = p, dpi=300, units="in", height=8, width=10)










