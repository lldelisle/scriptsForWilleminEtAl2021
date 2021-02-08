options(stringsAsFactors = F)
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
devtools::install_github("lldelisle/usefulLDfunctionsGR")
library(usefulLDfunctionsGR)
safelyLoadAPackageInCRANorBioconductor("GenomicRanges")
safelyLoadAPackageInCRANorBioconductor("seqinr")
# mm10 haploid and diploid
mm10.size <- c(2621487426, 5242974851)
names(mm10.size) <- c("haploid", "diploid")
sample.name <- "TgN3840_nCATS"
target.name <- "TgN3840_TLAderived"
size.target <- 64420
threshold.full.map.prop <- 0.99
threshold.full.map.diff <- 40
# Here we consider homozygous animals so it is like considering haploidy
ploidy.to.use <- "haploid"
mm10.size.to.use.for.enrichment <- mm10.size[ploidy.to.use]

# Get the length of the targetted reads from the fasta:
targetted.reads <- read.fasta("targetted_reads.fa", seqtype = "DNA")
seq.lengths <- sapply(targetted.reads, length)

# Get the mapping stats to evaluate the enrichment
map.stat <- read.delim(paste0(sample.name, "_mappingStats_det.txt"))
map.stat$sequenced <- map.stat$nMapped + map.stat$nUnmapped
map.stat$basesSeq <- map.stat$sequenced * map.stat$length
cat("The longest mapped read is:", max(map.stat$length[map.stat$nMapped > 0]), ".\n")
cat("The longest targetted read is:", max(seq.lengths), ".\n")

# Evaluate N50:
map.stat$cumBS <- cumsum(as.numeric(map.stat$basesSeq))
n50 <- map.stat$length[sum(map.stat$cumBS < sum(map.stat$basesSeq) / 2) + 1]
n50.targetted <- sort(seq.lengths)[sum(cumsum(sort(seq.lengths)) < sum(seq.lengths) / 2) + 1]
cat("The N50 is:", n50, ".\n")
cat("The N50 in targetted reads is:", n50.targetted, ".\n")


# Get the bed from the mapping on the construct
bed.file.gr <- grFromBedFile(paste0(sample.name, "_mapped_sorted_filtered2on", target.name, ".bed.gz"))
# Restrict to the interesting reads
bed.file.gr <- bed.file.gr[bed.file.gr$name %in% names(seq.lengths)]
bed.file.gr$map.length <- width(bed.file.gr)
bed.file.gr$read.length <- seq.lengths[bed.file.gr$name]
# Get the fully mapped reads:
full.map.reads.id <- with(mcols(bed.file.gr), which(map.length > read.length * threshold.full.map.prop | map.length >= read.length - threshold.full.map.diff))
# Check the one that support the insertion:
mapping.insertion <- as.data.frame(as.matrix(findOverlaps(bed.file.gr, GRanges(c('TgN3840_TLAderived'), IRanges(c(600, 64500), c(600, 64500))))))
mapping.insertion$map.length <- bed.file.gr$map.length[mapping.insertion$queryHits]
mapping.insertion$name <- bed.file.gr$name[mapping.insertion$queryHits]
mapping.insertion <- subset(mapping.insertion, name %in% mapping.insertion$name[duplicated(mapping.insertion$name)] )
mapping.insertion.ag <- aggregate(list(tot.map.length = mapping.insertion$map.length),
                                  by = list(name=mapping.insertion$name),
                                  FUN = sum)
mapping.insertion.ag$read.length <- seq.lengths[mapping.insertion.ag$name]
mapping.insertion.ag$map.length.min <- apply(mapping.insertion.ag[, c("read.length", "tot.map.length")], 1, min)
reads.insertion <- with(mapping.insertion.ag, name[which(tot.map.length > read.length * threshold.full.map.prop | tot.map.length >= read.length - threshold.full.map.diff)])

for (threshold in c(0, 2)){
  cat("Threshold:", threshold, "kb.\n")
  total.reads.threshold <- sum(map.stat$sequenced[map.stat$length > threshold * 1e3])
  cat("There are", length(seq.lengths[seq.lengths > threshold * 1e3]),
      "reads targetted out of", total.reads.threshold, ".\n")
  cat("There are", length(seq.lengths[seq.lengths > threshold * 1e3]),
      "reads targetted out of", sum(map.stat$nMapped[map.stat$length > threshold * 1e3]), " mapped.\n")
  fully.map.threshold <- names(full.map.reads.id[seq.lengths[names(full.map.reads.id)] > threshold * 1e3])
  insertion.threshold <- reads.insertion[seq.lengths[reads.insertion] > threshold * 1e3]
  cat("There are", length(fully.map.threshold),
      "reads targetted which fully mapped on the construct. And", length(insertion.threshold),
      "which support the insertion.\n")
  cat("Ratio fully mapped:", length(fully.map.threshold) / total.reads.threshold, ".\n")
  cat("Ratio both:", (length(fully.map.threshold) + length(insertion.threshold)) / total.reads.threshold, ".\n")
  nb.bases.fullymap.threshold <- sum(bed.file.gr$map.length[full.map.reads.id[fully.map.threshold]])
  nb.bases.insertion.threshold <- sum(mapping.insertion.ag$map.length.min[mapping.insertion.ag$name %in% insertion.threshold])
  cat("There are", format(nb.bases.fullymap.threshold, big.mark = ","),
      "bases from reads which fully mapped on the construct out of ", format(sum(map.stat$basesSeq[map.stat$length > threshold * 1e3]), big.mark = ","),
      ". And", format(nb.bases.insertion.threshold, big.mark = ","),
      "bases from reads which support the insertion.\n")
  cat("Ratio fully mapped:", nb.bases.fullymap.threshold / sum(map.stat$basesSeq[map.stat$length > threshold * 1e3]), ".\n")
  cat("Ratio both:", (nb.bases.fullymap.threshold + nb.bases.insertion.threshold) / sum(map.stat$basesSeq[map.stat$length > threshold * 1e3]), ".\n")
  cat("The construct is", format(size.target, big.mark = ","), "bp and a", ploidy.to.use, "genome is around", round(mm10.size.to.use.for.enrichment / 1e9, 1), "Gbp (",
      size.target / mm10.size.to.use.for.enrichment, ")\n")
  cat("So the enrichment is:", nb.bases.fullymap.threshold / sum(map.stat$basesSeq[map.stat$length > threshold * 1e3]) / (size.target / mm10.size.to.use.for.enrichment), "\n")
  cat("The construct with insertion is", format(size.target + 600, big.mark = ","), "bp and a", ploidy.to.use, "genome is around", round(mm10.size.to.use.for.enrichment / 1e9, 1), "Gbp (",
      (size.target + 600) / mm10.size.to.use.for.enrichment, ")\n")
  cat("So the enrichment is:", (nb.bases.fullymap.threshold + nb.bases.insertion.threshold) / sum(map.stat$basesSeq[map.stat$length > threshold * 1e3]) / ((size.target + 600) / mm10.size.to.use.for.enrichment), "\n")
  cat("\n\n")
}
