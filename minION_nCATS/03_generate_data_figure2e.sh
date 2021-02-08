#!/bin/bash

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 50G 
#SBATCH --cpus-per-task 1
#SBATCH --time 4:00:00
#SBATCH --job-name nCATSfig
#SBATCH --chdir /scratch/ldelisle/nCATS/

gitHubDirectory=$1

path="$PWD/"
sample=TgN3840_nCATS

pathResults=${path}/${sample}/
cd $pathResults

module purge
module load gcc/7.4.0
module load samtools/1.9
module load openblas/0.3.6-openmp
module load r/3.6.0

# Region of the transgene
samtools view ${sample}_mapped_sorted_filtered2.bam chr2:75123000-75160000 | cut -f 1 > targetted_reads.txt
# Region of insertion
samtools view ${sample}_mapped_sorted_filtered2.bam chr10:97018700-97019300 | cut -f 1 >> targetted_reads.txt
targetted_reads=$(cat targetted_reads.txt | tr "\n" ",")
zcat ${sample}_mapped_sorted_filtered2.bed.gz | awk -v ir=$targetted_reads 'BEGIN{
  split(ir,myreads,",")
  for (i in myreads){
    reads[myreads[i]] = 1
  }}
  {if($4 in reads){print}}' | gzip > ${sample}_mapped_sorted_filtered2_targetted.bed.gz

# samtools view -h $f | awk -v ir=$targetted_reads 'BEGIN{
#   split(ir,myreads,",")
#   for (i in myreads){
#     reads[myreads[i]] = 1
#   }}
#   {if($1~/^@/){print} else{if($1 in reads){print}}}' > ${sample}_mapped_sorted_filtered2_targetted.sam

# python ${gitHubDirectory}/scripts/minION/mappingStatsNanopore.py --output ${sample}_mappingStats_targetted.txt --outputDetail ${sample}_mappingStats_targetted_det.txt ${sample}_mapped_sorted_filtered2_targetted.sam

cat targetted_reads.txt | sort | uniq > targetted_reads_uniq.txt

# seqtk version 1.3.0
/home/ldelisle/softwares/seqtk/seqtk subseq combinedFastq/${sample}.fastq.gz targetted_reads_uniq.txt > targetted_reads.fastq
/home/ldelisle/softwares/seqtk/seqtk seq -A targetted_reads.fastq > targetted_reads.fa

# Generate stats:
Rscript ${gitHubDirectory}/scripts/minION/generateStats.R > ${gitHubDirectory}/minION_nCATS/supTableS4.txt

