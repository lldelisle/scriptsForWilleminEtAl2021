#!/bin/bash

#SBATCH -o slurm-%A_%2a.out
#SBATCH -e slurm-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 4G
#SBATCH --cpus-per-task 1
#SBATCH --time 2:00:00
#SBATCH --array=1
#SBATCH --job-name mergeBAM
#SBATCH --chdir /scratch/ldelisle/CNV

gitHubDirectory=$1

path="$PWD/"
pathForSraTable=${gitHubDirectory}/CopyNumber/sraTable.txt

module purge
module load gcc/7.4.0 #required for bowtie2, samtools, star and bedtools
module load samtools/1.9


sampleName=control_merged
samplesBamPath=$(cat $pathForSraTable | awk -v p="$path" '!($2~/TgN3840/){print p"/"$2"/"$2"_mapped_sorted_q30_rmdup.bam"}')
n=$(echo $samplesBamPath | wc -w)
sample="${sampleName}_neq$n"

pathResults=${path}/${sample}/
echo $sample
mkdir -p $pathResults
cd $pathResults

if [ ! -e ${sample}_mapped_sorted_q30_rmdup.bam ]; then
  samtools merge ${sample}_mapped_sorted_q30_rmdup.bam  $samplesBamPath
fi
if [ ! -e ${sample}_mapped_sorted_q30_rmdup.bam.bai ]; then
  samtools index ${sample}_mapped_sorted_q30_rmdup.bam
fi
if [ ! -e ${sample}_chr2_71_78.bam ]; then
  samtools view ${sample}_mapped_sorted_q30_rmdup.bam chr2:71000000-78000000 | awk -v OFS="\t" 'BEGIN{print "@HD\tVN:1.0\tSO:coordinate";print "@SQ\tSN:Fakechr2\tLN:7000000"}{$3="Fakechr2";$4=$4-71000000;print}' | samtools view -b - > ${sample}_chr2_71_78.bam
fi
