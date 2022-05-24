#!/usr/bin/env Rscript

library(ggplot2)
library(ggpubr)
library(wesanderson)


args = commandArgs(trailingOnly=TRUE)

d <- read.table(args[1], header=F, col.names=c("chrom", "start", "end", "ncov", "nhet", "snpcount"))
# Express heterozygosity in base pairs - filter out bins with less than 6000 well-covered sites
data <- subset(d, d$ncov >= 6000)
data$Heterozygosity <- round((data$snpcount/10000), digits = 2)

# Sort chromosomes
data$chrom <- as.numeric(as.factor(data$chrom))
sorted_data <- data[order(data$chrom),]


pal <- wes_palette("Zissou1", 100, type = "continuous")
# Plot genome-wide heterozygosity
#p_genome_wide_heterozygosity <- ggplot(data=sorted_data, aes(x=start, y=1))+facet_grid(chrom ~ ., switch='y')+geom_tile(aes(fill=Heterozygosity))+scale_fill_gradientn(colours = pal, breaks=seq(min(sorted_data$Heterozygosity),max(sorted_data$Heterozygosity),(max(sorted_data$Heterozygosity)-min(sorted_data$Heterozygosity))/4))+theme_void()+theme(text = element_text(size = 7))
p_genome_wide_heterozygosity <- ggplot(data=sorted_data, aes(x=start, y=1))+facet_grid(chrom ~ ., switch='y')+geom_tile(aes(fill=Heterozygosity))+scale_fill_gradientn(colours = pal, breaks=c(0.00,0.0150, 0.0300), limits=c(0.00, 0.0300))+theme_void()+theme(text = element_text(size = 7))


# Save figure
output_file <- sub('\\.heterozygosity.txt$', '', basename(args[1]))
fig_name <- paste(output_file, 'heterozygosity.pdf', sep = '.')
fig_name
ggsave(path = dirname(args[1]), filename = fig_name, p_genome_wide_heterozygosity, width = 18, height = 20, units = "cm")
