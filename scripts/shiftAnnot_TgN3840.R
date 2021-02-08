##This script takes as input an annotation file (bed, gtf...)
##It will shift the annotations.
##If an annotation overlap an insertion/deletion, it will be split in 2 annotations with the same names.
##Idem if it overlap an inversion.

options(scipen=999)
options(stringsAsFactors=F)
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
library(tools)

if(length(commandArgs(TRUE))==0){
  script.basename<-getSrcDirectory(function(x) {x})
  cat("Choose the file to convert.\n")
  fileToConvert<-file.choose()
  cat("Which is the number of the column with the chromosme name.\n")
  colChr<-as.numeric(readLines(con=stdin(),n=1))
  cat("Which is the number of the column with the start position.\n")
  colStart<-as.numeric(readLines(con=stdin(),n=1))
  cat("Which is the number of the column with the end position.\n")
  colEnd<-as.numeric(readLines(con=stdin(),n=1))
  cat("Which is the number of the column with the strand information (put 0 if there is no).\n")
  colStrand<-as.numeric(readLines(con=stdin(),n=1))
  cat("Full output path\n")
  outputName<-readLines(con=stdin(),n=1)
} else{
  if(commandArgs(TRUE)[1]=="-h" || commandArgs(TRUE)[1]=="--help"){
    cat("Usage: Rscript shiftAnnot_TgN3840.R fileToConvert colChr colStart colEnd [colStrand] fullOutputPath \n")
    stop()
  }
  fileToConvert<-commandArgs(TRUE)[1]
  colChr<-as.numeric(commandArgs(TRUE)[2])
  colStart<-as.numeric(commandArgs(TRUE)[3])
  colEnd<-as.numeric(commandArgs(TRUE)[4])
  if(length(commandArgs(TRUE)) > 5){
    colStrand<-as.numeric(commandArgs(TRUE)[5])
    outputName<-commandArgs(TRUE)[6]
  } else{
    colStrand<-0
    outputName<-commandArgs(TRUE)[5]
  }
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  script.basename <- dirname(script.name)
}
################################################################################
# Source the functions which are in the other file:
other.name <- file.path(script.basename, "shiftAnnotFunctions_compatibleInv.R")
source(other.name)

cat("Loading input file...")
annotationData <- usefulLDfunctions:::.readFileFromConditionOnNcols(fileToConvert,
                                                                    paste0(">=", max(colChr, colStart, colEnd, colStrand)),
                                                                    keepQuote=T)
cat("loaded.\n")

# Convert to UCSC format
if (!(grepl("chr", annotationData[1, colChr]))){
  annotationData[, colChr] <- paste0(rep("chr", nrow(annotationData)), annotationData[, colChr])
  annotationData[annotationData[, colChr] == "chrMT", colChr] <- "chrM"
  todelete <- grep("^chr(GL|JH)", annotationData[, colChr])
  annotationData <- annotationData[-todelete, ]
}

# First get data from the TgN38-40 fosmid

annot3840 <- getDFfromCoo(annotationData, colChr, colStart, colEnd, 
                        chrToGet="chr2", startToGet=75122684, endToGet=75160161)
if (nrow(annot3840) > 0){
  annot3840[, colChr] <- "chr10"
}

# And from CS39partial
annotCS39partial <- getDFfromCoo(annotationData, colChr, colStart, colEnd, 
                                 chrToGet="chr2", startToGet=75147232, endToGet=75147258)
if (nrow(annotCS39partial)>0){
  annotCS39partial[, colChr] <- "chr10"
}


# chr10 in TgN38-40 is:
# 97019824 bases WT
# 2034 bases of pEpiFos (97019825 to 97021858)
# TgN(38-40): 75122684 to 75160161 of chr2 wt (97021859 to 97059336)
# 7518 bp of pEpiFos not linearlized (97059337 to 97066854)
# TgN(38-40) partial (16143bp): 75122684 to 75138826 (97066855 to 97082997)
# CS39partial (27bp): 75147232 to 75147258 (97082998 to 97083024)
# 9bp from pEpiFos (97083025 to 97083033)
# bases from WT: 97019222 (from 97083034)

# Shift in chr2
brdf_chr2 <- read.delim(text = "genome\tbr1\tbr2\ttg1\tbr3\tbr4\ttg2
K4654\t75133816\t75133815\t95\t75133816\t75153815\t0")

shiftedAnnotations_chr2 <- shiftDFFromBR(annotationData, genome="K4654", brdf_chr2, colChr, colStart, colEnd,
                                         colStrand, verbose=T, chromoWithTg="chr2", splitIfOverlap=T)

# Construct chr10
annotations_chr10 <- annotationData[annotationData[, colChr] == "chr10", ]
shiftedAnnotations <- shiftedAnnotations_chr2[shiftedAnnotations_chr2[, colChr] != "chr10", ]

# First wt part
shiftedAnnotations <- rbind(shiftedAnnotations, 
                            getDFfromCoo(annotations_chr10, colChr, colStart, colEnd, 
                                         chrToGet="chr10", startToGet=1, endToGet=97019824))
# First complete TgN(38-40)
annot3840_shifted <- shiftDFOfBP(annot3840, colStart, colEnd, bpToShift=97021858)

# Second partial TgN(38-40)
annot3840_partial <- getDFfromCoo(annot3840, colChr, colStart, colEnd,
                                chrToGet="chr10", startToGet=1, endToGet=16143)
annot3840_partial_shifted <- shiftDFOfBP(annot3840_partial, colStart, colEnd, bpToShift=97066854)

# partial CS39
annotCS39partial_shifted <- shiftDFOfBP(annotCS39partial, colStart, colEnd, bpToShift=97082997)

# Last wt part
annot_end_wt <- getDFfromCoo(annotations_chr10, colChr, colStart, colEnd, 
                             chrToGet="chr10", startToGet=97019222, endToGet=1e15)

annot_end_wt_shifted <- shiftDFOfBP(annot_end_wt, colStart, colEnd, bpToShift=97083033)

shiftedAnnotations <- rbind(shiftedAnnotations, 
                            annot3840_shifted,
                            annot3840_partial_shifted,
                            annotCS39partial_shifted,
                            annot_end_wt_shifted)

dir.create(dirname(outputName), showWarnings=F)
write.table(shiftedAnnotations, file=outputName,
            sep='\t', row.names=F, col.names=F, quote=F)