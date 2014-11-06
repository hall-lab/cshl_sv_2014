#!/usr/bin/env Rscript

# Draw a histogram from a text file
args <- commandArgs(trailingOnly=TRUE)
file <- args[1]
x <- as.numeric(scan(file))

pdf(paste0(args[1], '.pdf'), height=4, width=5)
hist(x, col='steelblue3', breaks=50, main=paste0('Histogram of ', args[1]), xlab=args[1])
dev.off()