#!/bin/bash

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 10G
#SBATCH --cpus-per-task 1
#SBATCH --time 24:00:00
#SBATCH --job-name getFastq
#SBATCH --chdir /scratch/ldelisle/HiC

gitHubDirectory=$1

module purge
module load sra-toolkit/2.9.6

mkdir -p fastq
cd fastq
# Download fastq of this analysis
if [ ! -e ${gitHubDirectory}/HiC/sraTable.txt ]; then
  echo "sraTable.txt does not exists"
  exit 1
fi
while read line; do
  sra=$(echo ${line} | awk '{print $1}')
  output=$(echo ${line} | awk '{print $2}')
  fasterq-dump -o ${output}.fastq ${sra}
  mv ${output}_1.fastq ${output}_R1.fastq
  mv ${output}_2.fastq ${output}_R2.fastq
  gzip ${output}_R1.fastq
  gzip ${output}_R2.fastq
done < ${gitHubDirectory}/HiC/sraTable.txt

# Concatenate files from 2 round of resequencing
for f in *_hiseq*.fastq.gz; do
  cat ${f} ${f/hiseq/nextseq} > ${f/_hiseq/}
done

cd ..
# Download hicup0.7.3
if [ ! -e hicup_v0.7.3.tar.gz ]; then
  wget http://www.bioinformatics.babraham.ac.uk/projects/hicup/hicup_v0.7.3.tar.gz
fi
if [ ! -e hicup_v0.7.3 ]; then
  tar -zxvmf hicup_v0.7.3.tar.gz 
fi
# We adapt hicup v0.7.3 to make it work on the mutant genome
# which has duplicated sequences
if [ ! -e hicup_v0.7.3_mod ]; then
  cp -r hicup_v0.7.3 hicup_v0.7.3_mod
fi
if [ ! -e hicup_v0.7.3_mod/hicup_mapper.ori ]; then
  cp hicup_v0.7.3_mod/hicup_mapper hicup_v0.7.3_mod/hicup_mapper.ori
  sed -i "s/sub bowtie2MultiMapRead {/sub bowtie2MultiMapRead {\nreturn 0;/" hicup_v0.7.3_mod/hicup_mapper
fi
