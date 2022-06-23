#!/usr/bin/env Rscript

library(ggplot2)


args = commandArgs(trailingOnly=TRUE)

# Read in file with information on density of ancestral conserved elements in a certain N-bp window
data <- read.table(args[1], header = F, col.names=c("Chrom", "Start", "End", "Density", "Base_bairs", "Length", "Fraction"))
file_name <- basename(args[1])
directory <- dirname(args[1])

# Generate density plot, one for each chromosome
max <- max(data$Density)
chrOrder <-c(paste(1:80), "Z", "W")
data$Chrom <- factor(data$Chrom, chrOrder)
df <- data[order(data$Chrom), ]
plot <- ggplot(data=df, aes(x=Start, y=Density))+geom_histogram(stat="identity", color='red')+xlab('')+ylab('Density in 100-kb window')+facet_wrap(~Chrom, nrow=1, scales="free")+theme_bw()+theme(text=element_text(size=8), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+ylim(0, max)+ggtitle(args[2])
ggsave(path = directory, filename = paste(file_name, '.pdf', sep = ''), plot, width = 70, height = 7, units = "cm")


