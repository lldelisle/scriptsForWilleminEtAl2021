#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 40G 
#SBATCH --cpus-per-task 4
#SBATCH --time 20:00:00
#SBATCH --array=3-4
#SBATCH --job-name HiCUP
#SBATCH --chdir /scratch/ldelisle/HiC/

gitHubDirectory=$1

path="$PWD/"
pathForTableWithSamples="${gitHubDirectory}/HiC/table_HiC.txt"
pathForScripts=${gitHubDirectory}/scripts/
genome="mm10_TgN3840"
restSeq="^GATC"
restName="DpnII"
genomePath="/scratch/ldelisle/mutantGenome/TgN3840/"
pathForB2Index="${genomePath}/${genome}"
pathForFasta="${genomePath}/${genome}.fa"
pathForSizes="${genomePath}/${genome}.fa.fai"
pathForHiCUP="${path}/hicup_v0.7.3_mod/"
pathForFastq="${path}/fastq/"

bins="40"

nThreads=4

module purge
module load gcc/7.4.0  openblas/0.3.6-openmp
module load r/3.6.0
module load bowtie2/2.3.5
module load samtools/1.9
module load htslib/1.9

pathForR=$(which R)
pathForBowtie2=$(which bowtie2)


pathForDigest="${path}/${genome}_digester_${restName}.txt.gz"

#Digest 
if [ ! -e $pathForDigest ]; then
  ${pathForHiCUP}hicup_digester --re1 $restSeq,$restName --genome $genome --zip --outdir $path $pathForFasta
  mv ${path}/Digest_${genome}_${restName}* $pathForDigest
else
  # As it is run in an array, it is better to wait it is finished.
  sleep 5m
fi


sample=$(cat $pathForTableWithSamples | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}')
fastqFileR1=$(cat $pathForTableWithSamples | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $3}')
fastqFileR2=$(cat $pathForTableWithSamples | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $4}')

fullPathR1=${pathForFastq}/${fastqFileR1}
fullPathR2=${pathForFastq}/${fastqFileR2}


mkdir -p $sample
cd $sample


# Check if an output bam exists
inputBAM=$(find . -name "*.hicup.bam")

if [ -z $inputBAM ]; then
  # Run hicup
  ${pathForHiCUP}hicup --bowtie2 $pathForBowtie2 --digest $pathForDigest --format Sanger --index $pathForB2Index --keep --threads 28 --zip --r $pathForR $fullPathR1 $fullPathR2
  # Update the inputBAM variable
  inputBAM=$(find . -name "*.hicup.bam")
fi

pathForDigestNGZ="${path}/${genome}_digester_${restName}.txt"

if [ ! -e $pathForDigestNGZ ]; then
  gunzip -c $pathForDigest > $pathForDigestNGZ
fi

# Convert the bam to validPair in juicebox format
# The python script contrary to the provided converter
# Keep the fragment id and use the middle of the fragment
if [ ! -e ${sample}.validPairs_nonSorted.txt.gz ]; then
  python ${pathForScripts}/fromHicupToJuicebox.py  \
    --fragmentFile $pathForDigestNGZ --colForChr 1 --colForStart 2 \
    --colForEnd 3 --colForID 4 --lineToSkipInFragmentFile 2 \
    --useMid $inputBAM | gzip > ${sample}.validPairs_nonSorted.txt.gz
fi

inputValidPairs=${sample}.validPairs_nonSorted.txt.gz
pairs=${sample}.validPairs.csort.gz

if [ ! -e $pathForSizes ]; then
  samtools faidx $pathForFasta
fi

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
exists=$(conda info --envs | awk '$1=="hicexplorer3.6"{print}' | wc -l)
if [ $exists -ne 1 ]; then
  conda create -y -n hicexplorer3.6 hicexplorer=3.6 pygenometracks=3.6
fi
conda activate hicexplorer3.6
cooler --version
# cooler, version 0.8.10

if [ ! -e $pairs ]; then
  # sort and index the pairs with cooler and tabix
  cooler csort -i tabix -c1 3 -c2 7 -p1 4 -p2 8 -o $pairs $inputValidPairs $pathForSizes
fi

for bin in $bins; do
  if [ ! -e ${genome}.${bin}kb.bins ]; then
    cooler makebins $pathForSizes "${bin}000" > ${genome}.${bin}kb.bins
  fi
  if [ ! -e ${sample}.${bin}kb.cool ]; then
    cooler cload tabix -p 1 -c2 7 -p2 8 --assembly $genome ${genome}.${bin}kb.bins $pairs ${sample}.${bin}kb.cool
    cp ${sample}.${bin}kb.cool ${sample}_raw.${bin}kb.cool
    echo "Balancing"
    cooler balance --cis-only ${sample}.${bin}kb.cool
    echo "Balanced"
  fi
done
