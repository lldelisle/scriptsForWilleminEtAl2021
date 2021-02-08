#!/bin/bash

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 1G
#SBATCH --cpus-per-task 1
#SBATCH --time 1:00:00
#SBATCH --job-name pre4C

desc=$1
gitHubDirectory=$2
pathWithMutantGenome=$3

pathWithInstall="$PWD/"
pathWithTableWithGenomes="${gitHubDirectory}/mutantGenome/table.txt"
pathForScripts="${gitHubDirectory}/scripts/4C/"
pathWithBRFiles="${pathForScripts}/"
pathWithtemplate4cfile="${gitHubDirectory}/4C/template_4cfile.fa"


module purge
module load gcc/7.4.0
module load openblas/0.3.6-openmp
module load r/3.6.0

wd="${pathWithInstall}/${desc}"
if [ ! -e ${wd} ]; then
  echo "${wd} does not exists"
  exit 1
fi
cd ${wd}

# Count the number of genomes in the table
if [ ! -e ${pathWithTableWithGenomes} ]; then
  echo "${pathWithTableWithGenomes} does not exists"
  exit 1
fi
n=$(cat ${pathWithTableWithGenomes} | awk 'END{print NR}')
# For each genome
for ((i=1;i<=${n};i++)); do
  genome=$(cat ${pathWithTableWithGenomes} | awk -v i=${i} 'NR==i{print $1}')
  brFiles=$(cat ${pathWithTableWithGenomes} | awk -v i=${i} 'NR==i{print $2}')
  if [ -z ${brFiles} ]; then
    assembly=${genome}
    if [ ${genome} = "mm10" ]; then
      genome="Wt"
    fi
  else
    assembly=mm10_${genome}
  fi
  if [ -e $PWD/res_files_mapping_${genome} ]; then
    # If there are mapping data:
    echo "preparing ${assembly}"
    # Copy the template with a name which matches the genome:
    if [ ! -e ${pathWithtemplate4cfile} ]; then
      echo "${pathWithtemplate4cfile} does not exists"
      exit 1
    fi
    cp ${pathWithtemplate4cfile} template_4cfile_${genome}.fa
    if [ ! -z ${brFiles} ]; then
      # We need to shift the coordinates which were in the template
      if [ ! -e ${pathForScripts}/shiftPrimerFileWithMultipleBR_poolMultipleVP.sh ]; then
        echo "${pathForScripts}/shiftPrimerFileWithMultipleBR_poolMultipleVP.sh does not exists"
        exit 1
      fi      
      bash ${pathForScripts}/shiftPrimerFileWithMultipleBR_poolMultipleVP.sh ${brFiles} template_4cfile_${genome}.fa ${pathWithBRFiles} ${pathForScripts}
      # The new template_4cfile_${genome}.fa has the coordinates which match the genome
    fi
    # To run the 4c pipline the header of the 4c primer file needs to match each sample
    # We also need a configuration file which specifies how replicates are merged together and where are the bam files.
    if [ ! -e ${pathForScripts}/build_config_4Cseq_5Mb_and_primerfile_withReplicates.sh ]; then
      echo "${pathForScripts}/build_config_4Cseq_5Mb_and_primerfile_withReplicates.sh does not exists"
      exit 1
    fi    
    bash ${pathForScripts}/build_config_4Cseq_5Mb_and_primerfile_withReplicates.sh $PWD ${genome} ${assembly} ${pathWithMutantGenome} template_4cfile_${genome}.fa
  fi
done
