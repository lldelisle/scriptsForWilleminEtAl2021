#!/bin/bash -l

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 30G
#SBATCH --cpus-per-task 1
#SBATCH --time 24:00:00
#SBATCH --job-name preTLA
#SBATCH --chdir /scratch/ldelisle/TLA

# This script will create the bowtie2 index for the backbones
# Get the positions of the motif

gitHubDirectory=$1

path="$PWD/"
pathForIndexBackBones="${path}/backBones/"
pathForFastaBackBones="${gitHubDirectory}/TLA/fasta/"
pathForFasta="/home/ldelisle/genomes/fasta/"
genome=mm10
pathForScripts=${gitHubDirectory}/scripts/TLA/

module purge
module load gcc/7.4.0 #required for bowtie2, samtools
module load bowtie2/2.3.5
module load samtools/1.9
# Require fastx_toolkit
export PATH=$PATH:~/softwares/fastx_toolkit/bin/
# Require python with biopython
module load python/3.7.3


# backBone = DNA which contains sequences susceptible to have been integrated
# We assume they can be circular so we will duplicate them to not detect false positive
# At the circularized boundary
# If there are multiple fasta files, they can be space separated:
backBonesSpace=TgN3840_fosmid

# Build bowtie2 indices
mkdir -p $pathForIndexBackBones
for backBone in $backBonesSpace; do
  backBoneFasta=${pathForFastaBackBones}/${backBone}.fa
  fai=${backBoneFasta}.fai
  if [ ! -e ${fai} ]; then
    samtools faidx $backBoneFasta
  fi
  myTemplateBB="${pathForIndexBackBones}/${backBone}"
  if [ ! -e ${myTemplateBB}.fa ]; then
    cp $backBoneFasta ${myTemplateBB}.fa
  fi
  if [ ! -e ${myTemplateBB}.4.bt2 ]; then
    bowtie2-build ${myTemplateBB}.fa ${myTemplateBB}
  fi
done

# Also build bowtie2 indices for tail-to-head tail-to-tail head-to-head backbone
for backBone in $backBonesSpace; do
  backBoneFasta=${pathForFastaBackBones}/${backBone}_TH_TT_HH.fa
  if [ ! -e $backBoneFasta ]; then
    if [ -e ${pathForFastaBackBones}/${backBone}.fa ]; then
      # I build the rev complement
      cat ${pathForFastaBackBones}/${backBone}.fa | sed 1d | tr -d "\n" | rev | tr "ATCGatcg" "TAGCtagc" > ${backBone}_revcomp_temp.seq
      cat ${pathForFastaBackBones}/${backBone}.fa ${pathForFastaBackBones}/${backBone}.fa ${backBone}_revcomp_temp.seq ${pathForFastaBackBones}/${backBone}.fa | awk -v name=${backBone}_TH_TT_HH 'BEGIN{print ">"name}!($0~/>/){print}'> ${backBone}_TH_TT_HH_temp.fa
      # To be able to build a fasta index, we need to reformat this backbone_TH_TT_HH:
      fasta_formatter -i ${backBone}_TH_TT_HH_temp.fa -o $backBoneFasta -w 60
    else
      echo "${pathForFastaBackBones}/${backBone}.fa does not exists."
      exit 1
    fi
  fi
  fai=${backBoneFasta}.fai
  if [ ! -e ${fai} ]; then
    samtools faidx $backBoneFasta
  fi
  myTemplateBB="${pathForIndexBackBones}/${backBone}_TH_TT_HH"
  if [ ! -e ${myTemplateBB}.fa ]; then
    cp $backBoneFasta ${myTemplateBB}.fa
  fi
  if [ ! -e ${myTemplateBB}.4.bt2 ]; then
    bowtie2-build ${myTemplateBB}.fa ${myTemplateBB}
  fi
done

# Get RE positions
for backBone in $backBonesSpace $genome; do
  backBoneFasta=${pathForFastaBackBones}/${backBone}.fa
  what="${backBone}"
  if [ "$backBone" = "$genome" ]; then
    backBoneFasta=${pathForFasta}/${genome}.fa
    what=${genome}
  fi
  python ${pathForScripts}findAllREPos.py \
    --input ${backBoneFasta} \
    --output allCATG_${what}.txt --motif CATG --verbose
done
