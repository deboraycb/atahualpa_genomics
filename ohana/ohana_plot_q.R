#! /usr/bin/env Rscript

require(tidyr)
require(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

# file is the ohana q matrix output of qpas
file <- args[1]
# keyfile is a file without header, then one sample per line, first column is sample name, 2nd column is population name *in same order as ohana output*:
keyfile <- args[2]
# pop order is a file with the populations in the order you want them sorted in the plot
popfile <- args[3]

key <- read.table(keyfile, as.is=T, col.names=c('samples', 'pops'))
data <- read.table(file, skip=1)
data$sample <- key$samples
data$pop <- key$pops
poporder <- read.table(popfile, as.is=T)[,1]

data2 <- data
#to order pops in the plot the way I wanted (not necesary if the order doesn't matter to you):
data2$pop <- factor(data$pop, levels=poporder)
data2 <- data2[order(data2$pop),]
                    
data_long <- gather(data2, -sample, -pop, key = "Component", value = "Proportion")

print(aggregate(data2[,1:3], by=list(data2$pop), FUN=function(x) round(mean(x)*100, 1)))

colors <- c("V1"='#1f78b4',"V2"='#33a02c',"V3"='#e31a1c',"V4"='#ff7f00',"V5"='#6a3d9a',"V6"='#b15928',
    "V7"='#a6cee3',"V8"='#b2df8a')#,"V9"='#fb9a99',"V10"='#fdbf6f',"V11"='#cab2d6',"V12"='#ffff99',
#    "V13"='#4d6680',"V14"='#4d8051',"V15"='#804d57', "V16"='#604d80', "V17"='#806f4d')

plotname <- paste0(strsplit(file, "\\.")[[1]][1],
                           "_bar-plot_orig.svg")
ggplot(data_long) +   
    geom_col(aes(x=sample, y=Proportion, fill = Component)) +
    scale_fill_manual(values = colors)+
#                      guide=FALSE)+
    facet_grid(.~pop,scales = "free",space="free",
               switch="x") +
    theme(panel.spacing = unit(0, "lines"), 
          strip.background = element_rect(fill=NA, colour="gray"),
          axis.ticks = element_blank(), axis.text.x = element_blank(),
          strip.text.x = element_text(size=6)
          )+
    xlab("Population")
ggsave(plotname, h=2, w=6)
