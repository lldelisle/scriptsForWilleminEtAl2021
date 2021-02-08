#!/bin/bash

#SBATCH -o slurm-%A_%2a.out
#SBATCH -e slurm-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 30G
#SBATCH --array=1
#SBATCH --cpus-per-task 1
#SBATCH --time 24:00:00
#SBATCH --job-name mutantGenome
#SBATCH --chdir /scratch/ldelisle/mutantGenome/

gitHubDirectory=$1

path="$PWD/"
pathForScripts="${gitHubDirectory}/scripts/4C/"
# Here there is no br file
pathWithBRFiles="${pathForScripts}"
firstEnzymeName="Nla"
secondEnzymeName="Dpn"
firstEnzymeRE="CATG"
secondEnzymeRE="GATC"
length=30

module purge
module load gcc/7.4.0 #required for bowtie2, samtools, star and bedtools
module load bowtie2/2.3.5
module load samtools/1.9
module load bedtools2/2.27.1
module load openblas/0.3.6-openmp 
module load r/3.6.0
# Require fastx_toolkit
export PATH=$PATH:~/softwares/fastx_toolkit/bin/


genome=TgN3840
brFiles=../shiftAnnot_TgN3840.R

if [ -z ${brFiles} ]; then
  assembly=${genome}
  if [ ${genome} = "mm10" ]; then
    genome="Wt"
  fi
else
  assembly=mm10_${genome}
fi

mkdir -p ${genome}

cd ${genome}

if [ ! -e ${assembly}.fa ]; then
  # We get chr2 which is mutant:
  wget https://zenodo.org/record/3826913/files/chr2_delCS3840.fa.gz?download=1 -O chr2.fa.gz
  # We get chr10 which is also mutant:
  wget https://zenodo.org/record/4292337/files/chr10_TgN3840.fa.gz?download=1 -O chr10.fa.gz
  # We get other chromosomes from UCSC:
  for i in 1 {3..9} {11..19} X Y M; do
    wget "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/chr${i}.fa.gz" -O chr${i}.fa.gz
  done
  # We concatenate everything
  zcat chr*.fa.gz | fasta_formatter -w 60 > ${assembly}.fa
fi

if [ ! -e ${assembly}.fa.fai ]; then
  samtools faidx ${assembly}.fa
fi

if [ ! -e ${assembly}.rev.1.bt2 ]; then
  bowtie2-build ${assembly}.fa ${assembly} &
fi

if [ ! -e ${assembly}.json ]; then
  if [ ! -e ${pathForScripts}/addNewGenomeToGenRep.R ]; then
    echo "${pathForScripts}/addNewGenomeToGenRep.R does not exists"
    exit 1
  fi
  Rscript ${pathForScripts}/addNewGenomeToGenRep.R ${assembly} ${assembly}.fa.fai ./
fi

if [ ! -e ${assembly}_rmsk.bed ]; then
  if [ ! -e ${path}mm10_rmsk.bed.gz ]; then
    echo "There is no ${path}mm10_rmsk.bed.gz. This file is obtained through UCSC website in tools table browser variations and repeat."
    exit 1
  fi
  if [ -z ${brFiles} ]; then
    gunzip -c ${path}mm10_rmsk.bed.gz > ${assembly}_rmsk.bed
  else
    if [[ "$brFiles" = *".R" ]]; then
      Rscript ${pathWithBRFiles}${brFiles} ${path}mm10_rmsk.bed.gz 1 2 3 6 ${assembly}_rmsk.bed
    else
      brFilesSpace=$(echo ${brFiles} | tr "," " ")
      for f in ${brFilesSpace}; do
        if [ ! -e ${pathWithBRFiles}${f} ]; then
          echo "${pathWithBRFiles}${f} does not exists but is required to convert the rmsk."
          exit 1
        fi
        ln -s ${pathWithBRFiles}${f} .
      done
      if [ ! -e ${pathForScripts}/shiftBedWithMultipleBR.sh ]; then
        echo "${pathForScripts}/shiftBedWithMultipleBR.sh does not exists"
        exit 1
      fi
      bash ${pathForScripts}/shiftBedWithMultipleBR.sh ${brFiles} ${path}mm10_rmsk.bed.gz ${assembly}_rmsk.bed 1
    fi
  fi
fi

libraryName=library_${genome}_${firstEnzymeName}${secondEnzymeName}_${length}bps

if [ ! -e ${libraryName}_segmentInfos.bed ]; then
  bedRptMaskFile=${assembly}_rmsk.bed
  if [ ! -e ${bedRptMaskFile} ]; then
    echo "${bedRptMaskFile} does not exists"
    exit 1
  fi
  if [ ! -e ${pathForScripts}/getRestEnzymeOccAndSeq.pl ]; then
    echo "${pathForScripts}/getRestEnzymeOccAndSeq.pl does not exists"
    exit 1
  fi
  if [ ! -e ${pathForScripts}/manual_createLibrary.sh ]; then
    echo "${pathForScripts}/manual_createLibrary.sh does not exists"
    exit 1
  fi
  ln -s ${pathForScripts}/getRestEnzymeOccAndSeq.pl .
  cat ${assembly}.fa | awk '{if($0~/^>/){print}else{print toupper($0)}}' | bash ${pathForScripts}/manual_createLibrary.sh -i - -m ${firstEnzymeRE} -s ${secondEnzymeRE} -l ${length} -r ${bedRptMaskFile} -n ${libraryName} > ${libraryName}.log 2>${libraryName}.err
fi
wait
