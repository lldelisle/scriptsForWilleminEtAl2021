#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 30G
#SBATCH --cpus-per-task 1
#SBATCH --time 1:00:00
#SBATCH --array=1
#SBATCH --job-name FREEC
#SBATCH --chdir /scratch/ldelisle/CNV

gitHubDirectory=$1

path="$PWD/"
pathForTable="${gitHubDirectory}/CopyNumber/table.txt"

module purge
module load gcc/7.4.0 #required for bowtie2, samtools, star and bedtools
module load samtools/1.9
export PATH=$PATH:/home/ldelisle/softwares/FREEC-11.5/src/

sampleName=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}')
controlName=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}')

if [ ! -e  ${path}/${sampleName}/${sampleName}_chr2_71_78.bam ]; then
  if [ ! -e ${path}/${sampleName}/${sampleName}_mapped_sorted_q30_rmdup.bam.bai ]; then
    samtools index ${path}/${sampleName}/${sampleName}_mapped_sorted_q30_rmdup.bam
  fi
  samtools view ${path}/${sampleName}/${sampleName}_mapped_sorted_q30_rmdup.bam chr2:71000000-78000000 | awk -v OFS="\t" 'BEGIN{print "@HD\tVN:1.0\tSO:coordinate";print "@SQ\tSN:Fakechr2\tLN:7000000"}{$3="Fakechr2";$4=$4-71000000;print}' | samtools view -b - > ${path}/${sampleName}/${sampleName}_chr2_71_78.bam
fi

sampleBAM="${sampleName}_chr2_71_78.bam"
controlBAM="${controlName}_chr2_71_78.bam"

echo -e "Fakechr2\t7000000" > FakeChr2.size

mkdir -p ${sampleName}vs${controlName}_chr2_71_78
cd ${sampleName}vs${controlName}_chr2_71_78
ln -s ${path}/${sampleName}/$sampleBAM .
ln -s ${path}/${controlName}/$controlBAM .

for ws in 1 2; do
  ln -s $sampleBAM ${sampleBAM}_${ws}kb
  ln -s $controlBAM ${controlBAM}_${ws}kb

  echo "[general]" > config_${sampleName}vs${controlName}_${ws}kb.txt
  echo "chrLenFile = ../FakeChr2.size" >> config_${sampleName}vs${controlName}_${ws}kb.txt
  echo "ploidy = 2" >> config_${sampleName}vs${controlName}_${ws}kb.txt
  echo "breakPointThreshold = .8" >> config_${sampleName}vs${controlName}_${ws}kb.txt
  echo "window = ${ws}000" >> config_${sampleName}vs${controlName}_${ws}kb.txt
  echo "numberOfProcesses = 1" >> config_${sampleName}vs${controlName}_${ws}kb.txt
  echo "breakPointType = 4" >> config_${sampleName}vs${controlName}_${ws}kb.txt
  echo "BedGraphOutput=TRUE" >> config_${sampleName}vs${controlName}_${ws}kb.txt
  echo "[sample]" >> config_${sampleName}vs${controlName}_${ws}kb.txt
  echo "mateFile = ${sampleBAM}_${ws}kb" >> config_${sampleName}vs${controlName}_${ws}kb.txt
  echo "inputFormat = BAM" >> config_${sampleName}vs${controlName}_${ws}kb.txt
  echo "mateOrientation = 0" >> config_${sampleName}vs${controlName}_${ws}kb.txt
  echo "[control]" >> config_${sampleName}vs${controlName}_${ws}kb.txt
  echo "mateFile = ${controlBAM}_${ws}kb" >> config_${sampleName}vs${controlName}_${ws}kb.txt
  echo "inputFormat = BAM" >> config_${sampleName}vs${controlName}_${ws}kb.txt
  echo "mateOrientation = 0" >> config_${sampleName}vs${controlName}_${ws}kb.txt

  freec -conf config_${sampleName}vs${controlName}_${ws}kb.txt
done

for f in *.cpn; do
  echo $f
  cat $f | awk -v OFS="\t" 'NR==1{begin=$2;score=$3}NR>1{print "chr2\t"71000000+begin"\t"71000000+$2"\t"score;begin=$2;score=$3}' > ${f}_backToNormal.bedgraph
done

for f in *.BedGraph ; do
  echo $f
  cat $f | awk -v OFS="\t" 'NF==4{$1="chr2";$2=$2+71000000;$3=$3+71000000;print}' > ${f}_backToNormal.bedgraph
done

for f in *ratio.txt ; do
  echo $f
  cat $f | awk -v OFS="\t" 'NR==2{begin=$2;score=$4}NR>2{print "chr2\t"71000000+begin"\t"71000000+$2"\t"2^score;begin=$2;score=$4}' > ${f}_meanRatio_backToNormal.bedgraph
done

# For GEO renaming:
for ws in 1 2; do
  cp ${sampleName}_chr2_71_78.bam_${ws}kb_ratio.BedGraph_backToNormal.bedgraph ${sampleName}_${ws}kb_ratio.bedgraph
  cp ${sampleName}_chr2_71_78.bam_${ws}kb_ratio.txt_meanRatio_backToNormal.bedgraph ${sampleName}_${ws}kb_meanRatio.bedgraph
done
