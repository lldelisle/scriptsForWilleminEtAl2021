options(stringsAsFactors = F)

path.For.HiC <- commandArgs(TRUE)[1]
output <- commandArgs(TRUE)[2]

hicup.reports.fn <- list.files(path.For.HiC, recursive = T, pattern = "HiCUP.*report.*.txt", full.names = T)
hicup.reports <- as.data.frame(do.call(rbind, lapply(hicup.reports.fn, read.delim)))[, c(2, 18, 20, 28, 30)]
colnames(hicup.reports) <- gsub("_1$", "", colnames(hicup.reports))
hicup.reports$sample <- basename(dirname(hicup.reports.fn))
write.table(hicup.reports[, c("sample", setdiff(colnames(hicup.reports), "sample"))],
    file=output, quote=F, sep="\t", row.names = F)
cat("\n\nNote: The number of pairs is increased when mapped on TgN3840 as no filtering on mapping quality was performed.\n",
    file = output, append = T)
