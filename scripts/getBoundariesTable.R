options(stringsAsFactors = F)

library(rtracklayer)
library(plyr)

path.For.HiC <- commandArgs(TRUE)[1]
my.regions.bed <- commandArgs(TRUE)[2]
output <- commandArgs(TRUE)[3]

boundaries.files <- list.files(path.For.HiC, recursive=T, pattern=".40kb.*_boundaries.bed", full.names=T)
meta.df <- as.data.frame(do.call(rbind, strsplit(sapply(boundaries.files, basename), "\\.|kb")))[, c(1, 2, 4)]
colnames(meta.df) <- c("sample", "binsize", "windowsize")
meta.df$file <- boundaries.files
meta.df$sample.ws <- paste0(meta.df$sample, "_", meta.df$windowsize)
boundaries.gr <- lapply(boundaries.files , import)

my.regions.gr <- import(my.regions.bed)
my.regions.gr.extended <- flank(my.regions.gr, both=T, width=160000)

overlap.my.regions <- function(x){
    return(subsetByOverlaps(my.regions.gr.extended, x))
}
my.regions.gr.in.boundaries <- lapply(boundaries.gr, overlap.my.regions)

my.boundaries.df <- do.call(rbind, lapply(1:length(boundaries.files), function(i){
    data.frame(file=boundaries.files[i], name=my.regions.gr.in.boundaries[[i]]$name)
}))
my.boundaries.df.annotated <- merge(my.boundaries.df, meta.df)
write.table(as.matrix(table(my.boundaries.df.annotated$sample.ws, my.boundaries.df.annotated$name)),
    file=output, quote=F, sep="\t")
