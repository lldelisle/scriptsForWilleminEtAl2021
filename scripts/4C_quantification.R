options(stringsAsFactors = F)

library(rtracklayer)
library(plyr)
# Define functions

Read.Bigwig <- function(path_to_bw){
  bw <- as.data.frame(import.bw(path_to_bw))
  return(bw)
}

Read.Bedgraph.Or.Bed <- function(path_to_bdg_or_bed){
  bdg_or_bed <- read.table(path_to_bdg_or_bed, sep = "\t", fill = TRUE, header = FALSE, skip = 1, quote = "\"", stringsAsFactors = F) #of course, skip=1 is fine here, because only the first line of the bedgraph contains comments and therefore has to be ommited; but this is not necessarly always the case (Lucille created a more general function that keeps everything once 4 columns are detected)
  return(bdg_or_bed)
  #Condensed form: return(read.table(path_to_bdg_or_bed, sep="\t", fill = TRUE, header = FALSE, skip = 1, quote = "\"", stringsAsFactors = F))
}

Extract.Region.Of.Interest.From.Bed <- function(path_to_bed_of_regions, region_name){
  bed_of_region <- Read.Bedgraph.Or.Bed(path_to_bed_of_regions)
  region_of_interest <- subset(bed_of_region, bed_of_region[,4] == region_name)
  return(region_of_interest)
  #Condensed form: return(subset(Read.Bedgraph.Or.Bed(path_to_bed_of_regions), region_bed[,4] == region_name))
}

Sum.On.Interval.with.Overlap.Limit.Included.For.bw <- function(path_to_bw, path_to_bed_of_regions, region_name){
  region_of_interest <- Extract.Region.Of.Interest.From.Bed(path_to_bed_of_regions, region_name)
  bw <- Read.Bigwig(path_to_bw)
  data_to_quantify <- subset(bw, bw$seqnames == region_of_interest[1, 1] & bw$start <= region_of_interest[1, 3] & bw$end >= region_of_interest[1, 2]) #Everything that will overlap the limit will be included.
  result <- sum(data_to_quantify$score)
  return(result)
  #Condensed form: return(sum(subset(bedgraph, bedgraph[, 1] == region_of_interest[1, 1] & bedgraph[, 2] <= region_of_interest[1, 3] & bedgraph[, 3] >= region_of_interest[1, 2])$V4))
}

# Args:
file1.To.Quantify <- commandArgs(TRUE)[1]
file2.To.Quantify <- commandArgs(TRUE)[2]
bed.With.Regions <- commandArgs(TRUE)[3]
region1  <- commandArgs(TRUE)[4]
region2 <- commandArgs(TRUE)[5]

values <- data.frame(file = rep(c(file1.To.Quantify, file2.To.Quantify), 2),
                     region = rep(c(region1, region2), each=2))
values <- ddply(values, .(file, region), mutate,
                quantity = Sum.On.Interval.with.Overlap.Limit.Included.For.bw(file, bed.With.Regions, region))
values <- ddply(values, .(file), mutate,
                tot.quantity = sum(quantity))
values <- ddply(values, .(file, region), mutate,
                norm.quantity = quantity/tot.quantity)
values <- ddply(values, .(region), mutate,
                fc = (norm.quantity[file == file1.To.Quantify] - norm.quantity[file == file2.To.Quantify]) / norm.quantity[file == file2.To.Quantify])
print(unique(values[, c("region", "fc")]))
