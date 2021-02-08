#!/bin/bash

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 50G 
#SBATCH --cpus-per-task 10
#SBATCH --time 24:00:00
#SBATCH --job-name nCATS
#SBATCH --chdir /scratch/ldelisle/nCATS/

gitHubDirectory=$1

path="$PWD/"
sample=TgN3840_nCATS

nThreads=10

# First map on mm10:
genome="mm10"
pathForFasta="/home/ldelisle/genomes/fasta/${genome}.fa"
pathForMinimap2Index="/scratch/ldelisle/genomes/minimap2/${genome}.mmi"
pathForMinimap2IndexStorage="/work/updub/minimap2/${genome}.mmi"

mkdir -p $(dirname $pathForMinimap2Index)

module purge
module load gcc/7.4.0
module load samtools/1.9
module load bedtools2/2.27.1
module load python/3.7.3

pathResults=${path}/${sample}/

mkdir -p $pathResults
echo $sample
cd $pathResults

mkdir -p combinedFastq
fastq="${pathResults}/combinedFastq/${sample}.fastq.gz"
if [ -e $fastq ]; then
  echo "The fastq already exists"
else
  cat ${path}/guppy_output/*.fastq | gzip > $fastq
fi

if [ ! -e $pathForMinimap2Index ]; then
  if [ -e $pathForMinimap2IndexStorage ]; then
    cp -r $pathForMinimap2IndexStorage $pathForMinimap2Index
  else
    ~/softwares/minimap2-2.15_x64-linux/minimap2 -t $nThreads -d $pathForMinimap2Index $pathForFasta
    cp -r $pathForMinimap2Index $pathForMinimap2IndexStorage
  fi
fi

if [ ! -e ${sample}_algn.sam ]; then
  ~/softwares/minimap2-2.15_x64-linux/minimap2 -t $nThreads -ax map-ont $pathForMinimap2Index $fastq > ${sample}_algn.sam
fi

if [ ! -e ${sample}_mappingStats_det.txt ]; then
  python ${gitHubDirectory}/scripts/minION/mappingStatsNanopore.py --output ${sample}_mappingStats.txt --outputDetail ${sample}_mappingStats_det.txt ${sample}_algn.sam
fi

# Due to a design flaw, BAM does not work with CIGAR strings with >65535 operations (SAM and CRAM work). 
# So this may not work
# Here it worked
if [ ! -e ${sample}_mapped_sorted.bam ]; then
  samtools sort --threads $nThreads -o ${sample}_mapped_sorted.bam ${sample}_algn.sam
fi

if [ ! -e ${pathForFasta}.fai ]; then
  samtools faidx $pathForFasta
fi

# We remove the non primary alignments
if [ ! -e  ${sample}_mapped_sorted_filtered2.bam.bai ]; then
  samtools view -b ${sample}_mapped_sorted.bam --threads $nThreads -F0x104 > ${sample}_mapped_sorted_filtered2.bam
  samtools index ${sample}_mapped_sorted_filtered2.bam
fi

if [ ! -e ${sample}_mapped_sorted_filtered2.bed.gz ]; then
  bedtools bamtobed -i ${sample}_mapped_sorted_filtered2.bam | gzip > ${sample}_mapped_sorted_filtered2.bed.gz
fi

# Map on the construct built from TLA:
target=TgN3840_TLAderived
pathForFasta="${gitHubDirectory}/minION_nCATS/${target}.fa"
pathForMinimap2Index="${target}.mmi"

if [ ! -e $pathForMinimap2Index ]; then
  ~/softwares/minimap2-2.15_x64-linux/minimap2 -t $nThreads -d $pathForMinimap2Index $pathForFasta
fi

if [ ! -e ${sample}_algn_on${target}.sam ]; then
  ~/softwares/minimap2-2.15_x64-linux/minimap2 -t $nThreads -ax map-ont $pathForMinimap2Index $fastq > ${sample}_algn_on${target}.sam
fi

if [ ! -e ${sample}_mapped_sorted_on${target}.bam ]; then
  samtools sort --threads 16 -o ${sample}_mapped_sorted_on${target}.bam ${sample}_algn_on${target}.sam
  samtools index ${sample}_mapped_sorted_on${target}.bam
fi

if [ ! -e  ${sample}_mapped_sorted_filtered2on${target}.bam.bai ]; then
  samtools view -b ${sample}_mapped_sorted_on${target}.bam -F0x104 > ${sample}_mapped_sorted_filtered2on${target}.bam
  samtools index ${sample}_mapped_sorted_filtered2on${target}.bam
fi

if [ ! -e ${sample}_mapped_sorted_filtered2on${target}.bed.gz ]; then
  bedtools bamtobed -i ${sample}_mapped_sorted_filtered2on${target}.bam | gzip > ${sample}_mapped_sorted_filtered2on${target}.bed.gz
fi
