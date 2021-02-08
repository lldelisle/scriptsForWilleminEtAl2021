options(stringsAsFactors = F)

library(rtracklayer)
library(plyr)

path.For.HiC <- commandArgs(TRUE)[1]
output <- commandArgs(TRUE)[2]

boundaries.files <- list.files(path.For.HiC, recursive=T, pattern=".40kb.240kb_boundaries.gff", full.names=T)
samples <- sapply(boundaries.files, function(s){gsub(".40kb.240kb_boundaries.gff", "", basename(s))})
boundaries <- lapply(boundaries.files , readGFF)

boundaries.gr <- lapply(boundaries, makeGRangesFromDataFrame, keep.extra.columns=T)

my.regions.gr <- GRanges(c("chr2", "chr10"), IRanges(c(73640001,95480001), c(75800000,97880000)))

overlap.my.regions <- function(x){
    return(subsetByOverlaps(x, my.regions.gr))
}
my.boundaries.gr <- lapply(boundaries.gr, overlap.my.regions)

my.boundaries.df <- do.call(rbind, lapply(1:length(my.boundaries.gr), function(i){
    df <- as.data.frame(my.boundaries.gr[[i]])[, c("seqnames", "start", "end", "score", "delta", "pvalue")]
    df$sample <- samples[i]
    return(df)
}))
write.table(my.boundaries.df[order(my.boundaries.df$seqnames, my.boundaries.df$start, my.boundaries.df$sample), ],
    file=output, quote=F, sep="\t", row.names=F)
