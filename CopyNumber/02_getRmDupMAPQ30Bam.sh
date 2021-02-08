#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 30G
#SBATCH --cpus-per-task 16
#SBATCH --time 2:00:00
#SBATCH --array=3-9
#SBATCH --job-name map
#SBATCH --chdir /scratch/ldelisle/CNV

gitHubDirectory=$1

# This script remove adapters 
# make alignment with Bowtie2
# select MAPQ30 alignments
# Remove duplicates with Picard

path="$PWD/"
pathForFastq="${path}/fastq/"
pathForSraTable=${gitHubDirectory}/CopyNumber/sraTable.txt
pathForIndex="/home/ldelisle/genomes/bowtie2/"
pathForFasta="/home/ldelisle/genomes/fasta/"
genome=mm10

module purge
module load gcc/7.4.0 #required for bowtie2, samtools, star and bedtools
module load bowtie2/2.3.5
module load samtools/1.9
module load picard/2.19
# cutadapt is working with python 3.6.1 built with intel 17.0.2

indexPath=${pathForIndex}${genome}

sample=$(cat $pathForSraTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}')
fastqFile=${sample}.fastq.gz
adapterSeq=CTGTCTCTTATACACATCTCCGAGCCCACGAGAC

pathResults=${path}/${sample}/
mkdir -p $pathResults

echo $sample

cd $pathResults

mkdir -p ${path}/reports
if [ ! -e ${path}/reports/${sample}_report-cutadapt.txt ]; then
  if [ ! -e "${pathForFastq}${fastqFile}" ]; then
    echo "NO fastq"
    exit 1
  fi
  if [ -z $adapterSeq ]; then
    cutadapt -m 15 -j 16 -q 30 -o ${pathResults}${sample}-cutadapt.fastq.gz "${pathForFastq}${fastqFile}" > ${pathResults}${sample}_report-cutadapt.txt
  else
    cutadapt -m 15 -j 16 -a $adapterSeq -q 30 -o ${pathResults}${sample}-cutadapt.fastq.gz "${pathForFastq}${fastqFile}" > ${pathResults}${sample}_report-cutadapt.txt
  fi
  cp ${pathResults}${sample}_report-cutadapt.txt ${path}/reports/
fi

if [ ! -e ${path}/reports/${sample}_mapping_stats.txt ];then
  bowtie2 -p 16 -x $indexPath -U ${sample}-cutadapt.fastq.gz 2> ${pathResults}${sample}_mapping_stats.txt  | samtools view --threads 16 -Su - | samtools sort --threads 16 -o ${pathResults}${sample}_mapped_sorted.bam
  samtools index ${pathResults}${sample}_mapped_sorted.bam
  cp ${pathResults}${sample}_mapping_stats.txt ${path}/reports/
fi

if [ ! -e ${pathResults}${sample}_mapped_sorted_q30.bam ]; then
  samtools view --threads 16 -b ${pathResults}${sample}_mapped_sorted.bam -q 30 > ${pathResults}${sample}_mapped_sorted_q30.bam
fi

mkdir -p ${path}bam

if [ ! -e ${path}bam/${sample}_mapped_sorted_q30_rmdup.bam ]; then
  picard MarkDuplicates SORTING_COLLECTION_SIZE_RATIO=0.15 I=${pathResults}${sample}_mapped_sorted_q30.bam O=${pathResults}${sample}_mapped_sorted_q30_rmdup.bam M=${pathResults}${sample}_q30_rmdup.log REMOVE_DUPLICATES=true AS=true
#   picard  MarkDuplicates -SORTING_COLLECTION_SIZE_RATIO 0.15 -I ${pathResults}${sample}_mapped_sorted_q30.bam -O ${pathResults}${sample}_mapped_sorted_q30_rmdup.bam -M ${pathResults}${sample}_q30_rmdup.log -REMOVE_DUPLICATES true -AS true
  cp ${pathResults}${sample}_q30_rmdup.log ${path}/reports/
  cp ${pathResults}${sample}_mapped_sorted_q30_rmdup.bam ${path}bam/
fi
