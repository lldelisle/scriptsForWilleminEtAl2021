#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 30G
#SBATCH --cpus-per-task 5
#SBATCH --time 24:00:00
#SBATCH --array=1-2
#SBATCH --job-name TLA
#SBATCH --chdir /scratch/ldelisle/TLA

nThreads=5

gitHubDirectory=$1

# This script make alignment with Bowtie2
# Coverage with bedtools
# Find potential breakpoints using python scripts

path="$PWD/"
pathForFastq=${path}fastq/
pathForIndex="/home/ldelisle/genomes/bowtie2/"
pathForFasta="/home/ldelisle/genomes/fasta/"
pathForIndexBackBones="${path}/backBones/"
pathForFastaBackBones="${gitHubDirectory}/TLA/fasta/"
genome=mm10
pathForScripts=${gitHubDirectory}/scripts/TLA/
pathForSRATable=${gitHubDirectory}/TLA/sraTable.txt
# Satellite repeat masker can be obtained on
# https://genome.ucsc.edu/cgi-bin/hgTables
# mm10
# group: Variation and Repeats
# track: RepeatMasker
# filter: repClass match Satellite
# output format BED
pathForSatellite=${path}/mm10_rpmk_Satellite.bed.gz

if [ ! -e $pathForSatellite ]; then
  echo "satellite file does not exists"
  exit 1
fi

module purge
module load gcc/7.4.0 #required for bowtie2, samtools and bedtools
module load bowtie2/2.3.5
module load samtools/1.9
module load bedtools2/2.27.1
module load python/3.7.3
# cutadapt is version 1.16 using python 3.6 installed with pip

indexPath=${pathForIndex}${genome}
pathForChromSizes=${pathForFasta}${genome}.fa.fai

sample=$(cat $pathForSRATable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}')
fastqFile=${sample}.fastq.gz
# Put here comma separated backBones
backBones=TgN3840_fosmid
backBonesSpace=$(echo $backBones | tr "," " ")
# We assume the sample is blablabla_viewpoint with no _ in blablabla
viewpoint=$(echo $sample | awk -F '_' '{print $2}')

mkdir -p ${path}$sample

cd ${path}$sample

# First we remove Nextera adapters and bad quality reads
if [ ! -e ${sample}_cutadapt_report.txt ]; then
  cutadapt -j $nThreads -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -q 30 -m 25 -o ${sample}_cutadapt.fastq.gz $pathForFastq/$fastqFile > ${sample}_cutadapt_report.txt
fi
# Then we map end to end --very-sensitive and keep the unmapped
if [ ! -e ${sample}_mapped_sorted.bam ]; then  
  bowtie2 -p $nThreads -x ${indexPath} --very-sensitive -U ${sample}_cutadapt.fastq.gz --un-gz ${sample}_unaligned.fastq.gz 2> ${sample}_mapping_stats_endtoend.txt | samtools view --threads $nThreads -Su - | samtools sort --threads $nThreads -o ${sample}_mapped_sorted.bam
fi
# Filter for low mapping quality to avoid high coverage at repeat
if [ ! -e ${sample}_q30.bam ]; then
  samtools view --threads $nThreads -b ${sample}_mapped_sorted.bam -q 30 > ${sample}_q30.bam
fi
# Get the coverage of mapq30
if [ ! -e ${sample}_q30_cov.bw ]; then
  bedtools genomecov -ibam ${sample}_q30.bam -bg -trackline > ${sample}_q30_cov_UCSC.bdg
  if [ ! -e ${sample}_q30_cov_sorted.bedgraph ]; then
    cat ${sample}_q30_cov_UCSC.bdg | grep -v track | LC_COLLATE=C sort -k1,1 -k2,2n > ${sample}_q30_cov_sorted.bedgraph
  fi
  bedGraphToBigWig ${sample}_q30_cov_sorted.bedgraph ${pathForChromSizes} ${sample}_q30_cov.bw
fi
# Get the coverage per Mb of mapq30
if [ ! -e coverage_${sample}_q30_1Mb.txt ]; then
  samtools view ${sample}_q30.bam | awk '{print $3"\t"int($4/1000000)}' | uniq -c | awk '{print $2"\t"$3*1000000"\t"($3+1)*1000000"\t"$1}' > coverage_${sample}_q30_1Mb.txt
fi
# Sort the coverage to get higher coverage on top
if [ ! -e sorted_coverage_${sample}_q30_1Mb.txt ]; then
  cat coverage_${sample}_q30_1Mb.txt | sort -k4,4nr > sorted_coverage_${sample}_q30_1Mb.txt
fi

# Get the coverage on the backbones tail-to-head tail-to-tail head-to-head coverage:
for backBone in $backBonesSpace; do
  backBoneFasta=${pathForFastaBackBones}/${backBone}_TH_TT_HH.fa
  if [ ! -e $backBoneFasta ]; then
    echo "$backBoneFasta does not exists."
    exit 1
  fi
  fai=${backBoneFasta}.fai
  if [ ! -e ${fai} ]; then
    echo "$fai does not exists."
    exit 1
  fi
  myTemplateBB="${pathForIndexBackBones}/${backBone}_TH_TT_HH"
  if [ ! -e ${myTemplateBB}.4.bt2 ]; then
    echo "${myTemplateBB}.4.bt2 does not exists."
    exit 1
  fi
  if [ ! -e ${sample}_on${backBone}_TH_TT_HH_sorted.bam ]; then
    bowtie2 -p $nThreads -x ${myTemplateBB} --very-sensitive -U ${sample}_cutadapt.fastq.gz 2> ${sample}_on${backBone}_TH_TT_HH_mapping_stats.txt  | samtools view --threads $nThreads -Su - | samtools sort --threads $nThreads -o ${sample}_on${backBone}_TH_TT_HH_sorted.bam
    samtools index ${sample}_on${backBone}_TH_TT_HH_sorted.bam
  fi
  if [ ! -e ${sample}_on${backBone}_TH_TT_HH.bw ]; then
    bedtools genomecov -ibam ${sample}_on${backBone}_TH_TT_HH_sorted.bam -bg -trackline > ${sample}_on${backBone}_TH_TT_HH_cov_UCSC.bdg
    if [ ! -e ${sample}_on${backBone}_TH_TT_HH_cov_sorted.bedgraph ]; then
      cat ${sample}_on${backBone}_TH_TT_HH_cov_UCSC.bdg | grep -v track | LC_COLLATE=C sort --parallel=$nThreads -T . -k1,1 -k2,2n > ${sample}_on${backBone}_TH_TT_HH_cov_sorted.bedgraph
    fi
    bedGraphToBigWig ${sample}_on${backBone}_TH_TT_HH_cov_sorted.bedgraph $fai ${sample}_on${backBone}_TH_TT_HH.bw
  fi
  if [ ! -e ${sample}_on${backBone}_TH_TT_HH_summary.bed ]; then
    bedtools bamtobed -i ${sample}_on${backBone}_TH_TT_HH_sorted.bam | gzip > ${sample}_on${backBone}_TH_TT_HH.bed.gz
    zcat ${sample}_on${backBone}_TH_TT_HH.bed.gz |  sort -k1,1 -k2,2n -k3,3n | cut -f 1-3,6 | uniq -c | awk -v OFS="\t" '{print $2, $3, $4, "x"$1, $1, $5}' > ${sample}_on${backBone}_TH_TT_HH_summary.bed
  fi
done

# First approach:
currentFastq=${sample}_unaligned

for backBone in $backBonesSpace $genome; do
  myTemplateBB="${pathForIndexBackBones}/${backBone}"
  what="${backBone}"
  fai=${pathForFastaBackBones}/${backBone}.fa.fai
  if [ "$backBone" = "$genome" ]; then
    myTemplateBB=${indexPath}
    what=${genome}
    fai=${pathForChromSizes}
  fi
  # Do local alignment
  if [ ! -e ${currentFastq}_localMapped_on${what}_sorted.bam ]; then
    bowtie2 -p $nThreads -x ${myTemplateBB} --very-sensitive-local -U ${currentFastq}.fastq.gz 2> ${currentFastq}_localMapped_on${what}_mapping_stats.txt  | samtools view --threads $nThreads -Su - | samtools sort --threads $nThreads -o ${currentFastq}_localMapped_on${what}_sorted.bam
    samtools index ${currentFastq}_localMapped_on${what}_sorted.bam
  fi
  if [ ! -e ${currentFastq}_localMapped_on${what}.bw ]; then
    bedtools genomecov -ibam ${currentFastq}_localMapped_on${what}_sorted.bam -bg -trackline > ${currentFastq}_localMapped_on${what}_cov_UCSC.bdg
    if [ ! -e ${currentFastq}_localMapped_on${what}_cov_sorted.bedgraph ]; then
      cat ${currentFastq}_localMapped_on${what}_cov_UCSC.bdg | grep -v track | LC_COLLATE=C sort --parallel=$nThreads -T . -k1,1 -k2,2n > ${currentFastq}_localMapped_on${what}_cov_sorted.bedgraph
    fi
    bedGraphToBigWig ${currentFastq}_localMapped_on${what}_cov_sorted.bedgraph ${fai} ${currentFastq}_localMapped_on${what}.bw
  fi
  # Ignore reads containing CATG
  if [ ! -e ${currentFastq}_localMapped_on${what}_withoutCATG.bam ]; then
    samtools view -H ${currentFastq}_localMapped_on${what}_sorted.bam >  tmp.header
    samtools view ${currentFastq}_localMapped_on${what}_sorted.bam | awk '$10!~/CATG/{print}' | cat tmp.header - | samtools view  --threads 28 -b > ${currentFastq}_localMapped_on${what}_withoutCATG.bam
  fi
  if [ ! -e ${currentFastq}_localMapped_on${what}_withoutCATG.bw ]; then
    bedtools genomecov -ibam ${currentFastq}_localMapped_on${what}_withoutCATG.bam -bg -trackline > ${currentFastq}_localMapped_on${what}_withoutCATG_cov_UCSC.bdg
    if [ ! -e ${currentFastq}_localMapped_on${what}_withoutCATG_cov_sorted.bedgraph ]; then
      cat ${currentFastq}_localMapped_on${what}_withoutCATG_cov_UCSC.bdg | grep -v track | LC_COLLATE=C sort --parallel=$nThreads -T . -k1,1 -k2,2n > ${currentFastq}_localMapped_on${what}_withoutCATG_cov_sorted.bedgraph
    fi
    bedGraphToBigWig ${currentFastq}_localMapped_on${what}_withoutCATG_cov_sorted.bedgraph ${fai} ${currentFastq}_localMapped_on${what}_withoutCATG.bw
  fi
done

# Second approach:
# Split unmapped reads at CATG
if [ ! -e ${sample}_unaligned_CATG_split.fastq.gz ]; then
  perl ${pathForScripts}/splitOneFastqOnRe.pl --input ${sample}_unaligned.fastq.gz --output ${sample}_unaligned_CATG_split.fastq.gz --re_seq CATG
fi

# Keep only reads above 25bp
if [ ! -e ${sample}_unaligned_CATG_split_over25.fastq.gz ]; then
  cutadapt  -j $nThreads -q 30 -m 25 -o ${sample}_unaligned_CATG_split_over25.fastq.gz ${sample}_unaligned_CATG_split.fastq.gz > ${sample}_cutadapt2_report.txt
fi
# Map the split reads on end-to-end
if [ ! -e ${sample}_unaligned_CATG_split_mapped_sorted.bam ]; then  
  bowtie2 -p $nThreads -x ${indexPath} --very-sensitive -U ${sample}_unaligned_CATG_split_over25.fastq.gz --un-gz ${sample}_unaligned_CATG_split_no${genome}.fastq.gz 2> ${sample}_unaligned_CATG_split_mapping_stats.txt | samtools view --threads $nThreads -Su - | samtools sort  --threads $nThreads -o ${sample}_unaligned_CATG_split_mapped_sorted.bam
fi
# Merge both end-to-end bam files
if [ ! -e ${sample}_on${genome}_sorted.bam ]; then
  samtools merge ${sample}_on${genome}_sorted.bam ${sample}_mapped_sorted.bam ${sample}_unaligned_CATG_split_mapped_sorted.bam
  samtools index ${sample}_on${genome}_sorted.bam
fi
# Filter for low mapping quality to avoid high coverage at repeat
if [ ! -e ${sample}_on${genome}_sorted_q30.bam ]; then
  samtools view --threads $nThreads -b ${sample}_on${genome}_sorted.bam -q 30 > ${sample}_on${genome}_sorted_q30.bam
fi
# Get the coverage of mapq30
if [ ! -e ${sample}_on${genome}_q30_cov.bw ]; then
  bedtools genomecov -ibam ${sample}_on${genome}_sorted_q30.bam -bg -trackline > ${sample}_on${genome}_q30_cov_UCSC.bdg
  if [ ! -e ${sample}_on${genome}_q30_cov_sorted.bedgraph ]; then
    cat ${sample}_on${genome}_q30_cov_UCSC.bdg | grep -v track | LC_COLLATE=C sort -k1,1 -k2,2n > ${sample}_on${genome}_q30_cov_sorted.bedgraph
  fi
  bedGraphToBigWig ${sample}_on${genome}_q30_cov_sorted.bedgraph ${pathForChromSizes} ${sample}_on${genome}_q30_cov.bw
fi
# Get the coverage per Mb of mapq30
if [ ! -e coverage_${sample}_1Mb.txt ]; then
  samtools view ${sample}_on${genome}_sorted_q30.bam | awk '{print $3"\t"int($4/1000000)}' | uniq -c | awk '{print $2"\t"$3*1000000"\t"($3+1)*1000000"\t"$1}' > coverage_${sample}_1Mb.txt
fi
# Sort the coverage to get higher coverage on top
if [ ! -e sorted_coverage_${sample}_1Mb.txt ]; then
  cat coverage_${sample}_1Mb.txt | sort -k4,4nr > sorted_coverage_${sample}_1Mb.txt
fi
# Now for each backbone we will keep only the unmapped end-to-end reads
currentFastq=${sample}_unaligned_CATG_split_no${genome}
for backBone in $backBonesSpace; do
  backBoneFasta=${pathForFastaBackBones}/${backBone}.fa
  if [ ! -e $backBoneFasta ]; then
    echo "$backBoneFasta does not exists."
    exit 1
  fi
  fai=${backBoneFasta}.fai
  if [ ! -e ${fai} ]; then
    echo "$fai does not exists."
    exit 1
  fi
  myTemplateBB="${pathForIndexBackBones}/${backBone}"
  if [ ! -e ${myTemplateBB}.4.bt2 ]; then
    echo "${myTemplateBB}.4.bt2 does not exists."
    exit 1
  fi
  if [ ! -e ${currentFastq}no${backBone}.fastq.gz ]; then
    bowtie2 -p $nThreads -x ${myTemplateBB} --very-sensitive -U ${currentFastq}.fastq.gz --un-gz ${currentFastq}no${backBone}.fastq.gz 2> ${currentFastq}_on${backBone}_mapping_stats.txt  | samtools view --threads $nThreads -Su - | samtools sort --threads $nThreads -o ${currentFastq}_on${backBone}_sorted.bam
    samtools index ${currentFastq}_on${backBone}_sorted.bam
  fi
  currentFastq=${currentFastq}no${backBone}
done

# Now for each backBone + genome
# We will map locally
for backBone in $backBonesSpace $genome; do
  myTemplateBB="${pathForIndexBackBones}/${backBone}"
  what="${backBone}"
  fai=${pathForFastaBackBones}/${backBone}.fa.fai
  if [ "$backBone" = "$genome" ]; then
    myTemplateBB=${indexPath}
    what=${genome}
    fai=${pathForChromSizes}
  fi
  if [ ! -e ${currentFastq}_localMapped_on${what}_sorted.bam ]; then
    bowtie2 -p $nThreads -x ${myTemplateBB} --very-sensitive-local -U ${currentFastq}.fastq.gz 2> ${currentFastq}_localMapped_on${what}_mapping_stats.txt  | samtools view --threads $nThreads -Su - | samtools sort --threads $nThreads -o ${currentFastq}_localMapped_on${what}_sorted.bam
    samtools index ${currentFastq}_localMapped_on${what}_sorted.bam
  fi
  if [ ! -e ${currentFastq}_localMapped_on${what}.bw ]; then
    bedtools genomecov -ibam ${currentFastq}_localMapped_on${what}_sorted.bam -bg -trackline > ${currentFastq}_localMapped_on${what}_cov_UCSC.bdg
    if [ ! -e ${currentFastq}_localMapped_on${what}_cov_sorted.bedgraph ]; then
      cat ${currentFastq}_localMapped_on${what}_cov_UCSC.bdg | grep -v track | LC_COLLATE=C sort --parallel=$nThreads -T . -k1,1 -k2,2n > ${currentFastq}_localMapped_on${what}_cov_sorted.bedgraph
    fi
    bedGraphToBigWig ${currentFastq}_localMapped_on${what}_cov_sorted.bedgraph $fai ${currentFastq}_localMapped_on${what}.bw
  fi
done

# Write statistics
if [ ! -e ${sample}_summary.txt ]; then
  n=$(cat ${sample}_cutadapt_report.txt | grep "Total reads processed" | awk '{print $NF}')
  echo -e "Input reads\t$n" > ${sample}_summary.txt
  n=$(cat ${sample}_cutadapt_report.txt | grep "Reads written" | awk '{print $(NF - 1)}')
  echo -e "Reads without adapters\t$n" >> ${sample}_summary.txt
  n=$(cat ${sample}_mapping_stats_endtoend.txt | grep "aligned 0 times" | awk '{print $1}')
  echo -e "Reads do not map end-to-end on ${genome}\t$n" >> ${sample}_summary.txt
  n=$(cat ${sample}_cutadapt2_report.txt | grep "Total reads processed" | awk '{print $NF}')
  echo -e "Reads obtained after split at CATG\t$n" >> ${sample}_summary.txt
  n=$(cat ${sample}_cutadapt2_report.txt | grep "Reads written" | awk '{print $(NF - 1)}')
  echo -e "Reads longer than 25bp\t$n" >> ${sample}_summary.txt
  n=$(cat ${sample}_unaligned_CATG_split_mapping_stats.txt | grep "aligned 0 times" | awk '{print $1}')
  echo -e "Reads do not map end-to-end on ${genome}\t$n" >> ${sample}_summary.txt
  currentFastq=${sample}_unaligned_CATG_split_no${genome}
  for backBone in $backBonesSpace; do
    n=$(cat ${currentFastq}_on${backBone}_mapping_stats.txt | grep "aligned 0 times" | awk '{print $1}')
    echo -e "Reads do not map end-to-end on ${backBone}\t$n" >> ${sample}_summary.txt
    currentFastq=${currentFastq}no${backBone}
  done
  for backBone in $backBonesSpace $genome; do
    what="${backBone}"
    if [ "$backBone" = "$genome" ]; then
      what=$genome
    fi
    n=$(cat ${currentFastq}_localMapped_on${what}_mapping_stats.txt | awk 'END{print $1}')
    echo -e "Mapping rate with local on $what\t$n" >> ${sample}_summary.txt
  done
fi

# Now look for potential break points
for backBone in $backBonesSpace $genome; do
  what="${backBone}"
  if [ "$backBone" = "$genome" ]; then
    what=$genome
  fi
  fullFile=${currentFastq}_findBP_full_on${what}.txt
  # First identify all potential Break Points (positions where reads aligned only partially from or to)
  if [ ! -s $fullFile ]; then
    python ${pathForScripts}findBP.py --input ${currentFastq}_localMapped_on${what}_sorted.bam --output ${currentFastq}_findBP_on${what}.txt --outputFull $fullFile --minCov 2 --minCovTot 5
  fi
  # Sort them by decreasing coverage
  if [ ! -s 01_sorted_$fullFile ]; then
    sort -k4,4nr $fullFile > 01_sorted_$fullFile
  fi
  # Filter for redundancy:
  # Breakpoints which differ from only one base of a more frequently found breakpoint.
  if [ ! -s 02_noRed_${fullFile} ]; then
    python ${pathForScripts}filter_forRedundancy.py --input 01_sorted_${fullFile} --output 02_noRed_${fullFile}
  fi
  # Filter for breakpoints at satellites or polyA or very close to viewpoint.
  if [ ! -s 03_rmBedACo_${fullFile} ]; then
    if [ "$backBone" = "${genome}" ]; then
      otherBed=${pathForSatellite}
    else
      otherBed=""
    fi
    if [ -n "$viewpoint" ]; then
      python ${pathForScripts}filter_forBedAndCo.py --input 02_noRed_${fullFile} --output 03_rmBedACo_${fullFile} --polyA --bed ${gitHubDirectory}/TLA/${viewpoint}.bed $otherBed
    else
      echo "THERE IS NO VIEWPOINT"
      if [ -n "$otherBed" ]; then
        python ${pathForScripts}filter_forBedAndCo.py --input 02_noRed_${fullFile} --output 03_rmBedACo_${fullFile} --polyA --bed $otherBed
      else
        python ${pathForScripts}filter_forBedAndCo.py --input 02_noRed_${fullFile} --output 03_rmBedACo_${fullFile} --polyA
      fi
    fi
  fi
  # Filter for reads where there was a sequencing error at the
  # level of the motif which prevented a good split
  if [ ! -s 04_rmMotif_${fullFile} ]; then
    python ${pathForScripts}filter_forMotifError.py --input 03_rmBedACo_${fullFile} --output 04_rmMotif_${fullFile} --motifFile ${path}/allCATG_${what}.txt
  fi
  # Write the consensus sequence of the breakpoint
  if [ ! -s consensusOn04_${fullFile} ]; then
    python ${pathForScripts}consensus_bp_from_full.py --input 04_rmMotif_${fullFile} --minCoverage 3 --output consensusOn04_${fullFile}
  fi
  # Get all reads well aligned on the breakpoint
  if [ ! -s readsAlignedOn04_${fullFile} ]; then
    python ${pathForScripts}check_bp_from_full.py --input 04_rmMotif_${fullFile} --output readsAlignedOn04_${fullFile}
  fi
done
