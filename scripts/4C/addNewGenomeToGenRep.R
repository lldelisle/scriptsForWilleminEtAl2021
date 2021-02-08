options(stringsAsFactors=F)
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("jsonlite")
library(tools)
################################################################################
if(length(commandArgs(TRUE))>0){
  if(commandArgs(TRUE)[1]=="-h" || commandArgs(TRUE)[1]=="--help"){
    cat("Usage: Rscript addNewGenomeToGenRepSimplifiedWithArgs.R assemblyName pathForChrSize outputFolder\n")
    stop()
  }
  assembly<-commandArgs(TRUE)[1]
  chromosomeSizeFile<-commandArgs(TRUE)[2]
  outputFolder<-commandArgs(TRUE)[3]
} else{
  connection<-stdin()
  cat("What is the name of your assembly: \n")
  assembly<-readLines(con=connection,1)
  cat("Select the file with chromosome size information.\n")
  pathForBr<-file.choose()
  cat("Where do you want to have the results (put the absolute path): \n")
  outputFolder<-readLines(con=connection,1)
}
################################################################################

newLine<-data.frame(assembly_type_id =NA, assembly_type_name =NA, bbcf_valid =NA,
                    created_at =format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ"), genome_id =1,
                    gtf_convention =NA, id =1, md5 =assembly, name =assembly,
                    nr_assembly_id =1, source_id =NA, source_name =NA,
                    ucsc_view_url =NA, updated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ"))
#You need a json with only your assembly
#To keep the structure, you need to create a dataframe with only one column which is a dataframe.
#Here is the trick:
ddfNL<-list(assembly=newLine)
class(ddfNL)<-"data.frame"
rownames(ddfNL)<-seq_len(nrow(newLine))
dir.create(outputFolder,showWarnings = F,recursive = T)
cat(toJSON(ddfNL,na="null"),file = paste0(outputFolder,"/",assembly,".json"))
##################
####Chromosomes###
##################
sizes<-read.delim(chromosomeSizeFile,h=F)
colnames(sizes)<-c("original_chr_name","length")
# You need to give an id to each chromosome
sizes$id<-sizes$original_chr_name
# The name does not have chr, it will be removed if it is UCSC format.
sizes$name<-gsub("^chr","",sizes$original_chr_name)
# The num is numeric so X, Y etc. should be changed.
# I decide to give the chromsomes with number the num corresponding
sizes$num<-sapply(sizes$name,function(n){tryCatch(as.numeric(n),warning=function(w){NA})})
sizes<-sizes[order(sizes$num),]
# Then I put the letters
lettersChr<-which(sizes$name%in%LETTERS)
lastNumberUsed <- tryCatch(max(sizes$num,na.rm = T), warning= function (w){0})
sizes$num[lettersChr[order(sizes$name[lettersChr])]]<-seq(from=lastNumberUsed+1,length.out = length(lettersChr))
# If there are additionnal contigs
# They are sorted by the chromosome they belong to
sizes$chromOri<-sapply(sizes$name,function(n){tryCatch(as.numeric(strsplit(n,"_")[[1]][1]),warning=function(w){NA})})
sizes<-sizes[order(sizes$num,sizes$chromOri,sizes$name),]
lastNumberUsed <- tryCatch(max(sizes$num,na.rm = T), warning= function (w){0})
sizes$num[is.na(sizes$num)]<-lastNumberUsed+(1:(sum(is.na(sizes$num))))
temp.chr_names<-apply(sizes,1,function(v){smallDF<-data.frame(assembly_id=1,chromosome_id=v["id"],
                                                              updated_at=format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ"),value=v["original_chr_name"])
                                          rownames(smallDF)<-NULL
                                          dsmallDF<-list(chr_name=smallDF)
                                          class(dsmallDF)<-"data.frame"
                                          rownames(dsmallDF)<-1
                                          return(dsmallDF)})
# Again the trick
df.temp.df1<-list(chr_names=temp.chr_names)
class(df.temp.df1)<-"data.frame"
rownames(df.temp.df1)<-seq_len(nrow(sizes))
# Then I add all the other fields
df.temp.df1$chr_type_id<-1
df.temp.df1$circular<-F
df.temp.df1$created_at<-format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ")
df.temp.df1$genome_id<-1
df.temp.df1$gi_number<-NA
df.temp.df1$id<-sizes$id
df.temp.df1$length<-sizes$length
df.temp.df1$md5<-NA
df.temp.df1$name<-sizes$name
df.temp.df1$num<-sizes$num
df.temp.df1$public<-T
df.temp.df1$refseq_locus<-NA
df.temp.df1$refseq_version<-NA
df.temp.df1$synonyms<-""
df.temp.df1$updated_at<-format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ")
# A json is generated with chromosomes of the assembly
ddfCA<-list(chromosome=df.temp.df1)
class(ddfCA)<-"data.frame"
rownames(ddfCA)<-seq_len(nrow(df.temp.df1))
cat(toJSON(ddfCA,na="null"),file = paste0(outputFolder,"/",assembly,"_chromosomes.json"))
# For the demultiplexing step, you need to put the 2 files in nr_assemblies/info_json
# For the mapping you also need to put your bowtie2 indexes in nr_assemblies/bowtie2 with the name of your assembly for example in my case galGal6.1.bt2 etc...
# For the 4C seq you need to generate a library. Do not forget to put your fasta in capital letter before.
