#!/bin/bash

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 50G
#SBATCH --cpus-per-task 1
#SBATCH --time 6:00:00
#SBATCH --job-name TADs_diff
#SBATCH --chdir /scratch/ldelisle/HiC

gitHubDirectory=$1

path="$PWD/"
bins="20 40"
sizes="240 320 480 800"
pathForTableWithSamples="${gitHubDirectory}/HiC/table_HiC.txt"

nbOfThreads=1

module purge

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
exists=$(conda info --envs | awk '$1=="hicexplorer3.6"{print}' | wc -l)
if [ $exists -ne 1 ]; then
  conda create -y -n hicexplorer3.6 hicexplorer=3.6 pygenometracks=3.6
fi
conda activate hicexplorer3.6

mkdir -p TADs_diff

cd TADs_diff

samples=$(cat $pathForTableWithSamples | awk '{print $2}')

for bin in $bins; do
  for size in $sizes; do
    for sample in ${samples}; do
      coolFile=${path}/${sample}/${sample}.${bin}kb.cool
      if [ ! -e ${sample}.${bin}kb.${size}kb_domains.bed ]; then
        echo "finding tads for $sample $bin"
        hicFindTADs -m ${coolFile} --minBoundaryDistance 200000 --correctForMultipleTesting bonferroni \
          --outPrefix ${sample}.${bin}kb.${size}kb \
          --minDepth ${size}000 --maxDepth $((2 * ${size}000)) --step $((2 * ${size}000))
      fi
    done
  done
done

# We will perform the difference between
mutant=E12_Limbs_TgN3840_map_TgN3840
control=E12_Limbs_Wt_map_TgN3840
for bin in $bins; do
  mutant_cool=${path}/${mutant}/${mutant}.${bin}kb.cool
  control_cool=${path}/${control}/${control}.${bin}kb.cool
  if [ ! -e E12_Limbs_TgN3840MinusWt_map_TgN3840.${bin}kb.cool ]; then
    echo "building E12_Limbs_TgN3840MinusWt_map_TgN3840.${bin}kb.cool"
    hicCompareMatrices -m ${mutant_cool} ${control_cool} \
      --outFileName E12_Limbs_TgN3840MinusWt_map_TgN3840.${bin}kb.cool \
      --operation diff
  fi
done
