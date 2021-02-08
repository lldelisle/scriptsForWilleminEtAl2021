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
#SBATCH --chdir /scratch/ldelisle/CNV

gitHubDirectory=$1

module purge
module load sra-toolkit/2.9.6

mkdir -p fastq
cd fastq
# Download fastq of this analysis
if [ ! -e ${gitHubDirectory}/CopyNumber/sraTable.txt ]; then
  echo "sraTable.txt does not exists"
  exit 1
fi
while read line; do
  sra=$(echo ${line} | awk '{print $1}')
  output=$(echo ${line} | awk '{print $2".fastq"}')
  fasterq-dump -o ${output} ${sra}
  gzip ${output}
done < ${gitHubDirectory}/CopyNumber/sraTable.txt
