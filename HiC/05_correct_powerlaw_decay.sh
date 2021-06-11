#!/bin/bash

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 50G
#SBATCH --cpus-per-task 1
#SBATCH --time 6:00:00
#SBATCH --job-name powerlaw_decay
#SBATCH --chdir /scratch/ldelisle/HiC

gitHubDirectory=$1

path="$PWD/"

nbOfThreads=1

module purge
module load gcc/7.4.0
module load openblas/0.3.6-openmp 
module load r/3.6.0
module load samtools/1.9

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
exists=$(conda info --envs | awk '$1=="hicexplorer_dev"{print}' | wc -l)
if [ $exists -ne 1 ]; then
  conda env create -f ${gitHubDirectory}/HiC/hicexplorer_dev.yml
fi
conda activate hicexplorer_dev

mkdir -p decay
cd decay

# Integration borders:
echo -e "chr10\t97019221\t97019222\trightBP\nchr10\t97019824\t97019825\tleftBP" > integration_borders.bed
Rscript ${gitHubDirectory}/scripts/shiftAnnot_TgN3840.R integration_borders.bed 1 2 3 integration_borders_TgN3840_temp.bed
cat integration_borders_TgN3840_temp.bed | awk '$4=="leftBP"{if($2==$3){print $1"\t"$2"\t"$3+1"\t"$4}else{print}}' > integration_borders_TgN3840.bed

sample=E12_Limbs_Wt_map_mm10
cool_file=${path}/${sample}/${sample}.40kb.cool

# First we adjust TADs boundaries so they correspond to fix bins:
domains=${path}/TADs_diff/${sample}.40kb.240kb_domains.bed
cat $domains | awk -v OFS="\t" -v bs=40000 '$2%bs!=0{$2 = int($2/bs) * bs}$3%bs!=0{$3 = (int($2/bs) + 1)* bs}{print}' > $(basename $domains .bed)_adjusted.bed
domains=$(basename $domains .bed)_adjusted.bed
# Get the coordinates of the Btg1 TAD
bedtools intersect -a $domains -b integration_borders.bed -wa -u > welcoming_tad.bed
# Compute average interaction frequency per genomic distance
# Genome wide
hicPlotDistVsCounts -m ${cool_file} --outFileData decay_genome_wide.txt -o decay_genome_wide.png
# On all TADs
hicPlotDistVsCounts -m ${cool_file} --outFileData decay_tads.txt -o decay_tads.png \
  --domains ${domains}
# On Btg1 TAD
hicPlotDistVsCounts -m ${cool_file} --outFileData decay_welcoming_tad.txt -o decay_welcoming_tad.png \
  --domains welcoming_tad.bed
# On the first TADs:
# Note: 4, 5, were excluded because they have NA bins
for i in 1 2 3 6 7 8; do
  awk -v i=$i 'NR==i{print}' $domains > tad_${i}.bed
  if [ ! -e decay_tad_${i}.txt ]; then
    hicPlotDistVsCounts -m ${cool_file} --outFileData decay_tad_${i}.txt \
      -o decay_tad_${i}.png \
      --domains tad_${i}.bed
  fi
done

Rscript $gitHubDirectory/scripts/hic_decay_all_weight.R FigS5A.pdf

# Now we apply the correction:
# On WT map on mutant
sample=E12_Limbs_Wt_map_TgN3840
# We get the coordinate of the TAD on mutant genome
domains=${path}/TADs_diff/${sample}.40kb.240kb_domains.bed
bedtools intersect -a $domains -b integration_borders_TgN3840.bed -wa -u > welcoming_tad_TgN3840.bed
welcoming_tad_start=$(cat welcoming_tad_TgN3840.bed | awk '{print $2}')
welcoming_tad_end=$(cat welcoming_tad_TgN3840.bed | awk '{print $3}')
plotting_end=99880000
# Get the alpha_value from the R analysis
alpha_value=$(cat powerlaw_decay_on_welcoming_weight.txt)
# Create an artificial chr10 to have smaller files:
echo -e "chr10\t${plotting_end}" > small_chr10.size
# Correct and store in new cool file
cooler dump -b --join ${path}/${sample}/${sample}.40kb.cool \
  -r chr10:0-${plotting_end} | awk -v OFS="\t" -v alpha=$alpha_value \
  -v startI=$(awk 'NR==1{print $2}' integration_borders_TgN3840.bed) \
  -v endI=$(awk 'NR==2{print $2}' integration_borders_TgN3840.bed) \
  -v startT=$welcoming_tad_start \
  -v endT=$welcoming_tad_end '
$2 < startI && $2 >= startT && $5 > endI && $5 < endT {
  # Correct of powerlaw decay:
  $8 *= (($5 - $2) / ($5 - $2 - 63812)) ** alpha
}
{
  print $1,$2,$3,$4,$5,$6,$8
}' | \
  cooler load --count-as-float --format bg2 \
      "small_chr10.size:40000" - ${sample}_welcoming_tad.cool

# For the mutant, just extract:
sample=E12_Limbs_TgN3840_map_TgN3840
cooler dump -b --join ${path}/${sample}/${sample}.40kb.cool \
  -r chr10:0-${plotting_end} | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$8}' | \
  cooler load --count-as-float --format bg2 \
      "small_chr10.size:40000" - ${sample}_welcoming_tad.cool


# Compute the difference bin per bin
# Outside the insertion region
mutant=E12_Limbs_TgN3840_map_TgN3840
control=E12_Limbs_Wt_map_TgN3840
awk '
NR==FNR{
  # Control:
  diff[$1"_"$2] -= $3
}
NR!=FNR{
  # Mutant:
  if (diff[$1"_"$2] != 0){
    # We only add the mutant value when wt is not 0
    # Which corresponds to the insertion region
    diff[$1"_"$2] += $3
  }
}
END{
  for(bin1bin2 in diff){
    split(bin1bin2, a, "_")
    print a[1]"\t"a[2]"\t"diff[bin1bin2]
  }
}' <(cooler dump ${control}_welcoming_tad.cool) \
  <(cooler dump ${mutant}_welcoming_tad.cool) | \
  cooler load --count-as-float --format coo \
      "small_chr10.size:40000" - ${mutant}Minus${control}_welcoming_tad.cool
