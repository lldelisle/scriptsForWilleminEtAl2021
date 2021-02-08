#! /bin/bash

# Define all paths
gitHubDirectory="/home/ldelisle/softwares/scriptsForWilleminEtAl2021/"
pathForBBCFutils="/home/ldelisle/softwares/" # BBCFutils is found at git clone https://github.com/bbcf/bbcfutils.git -b standalone
mutantGenomeDirectory="/scratch/ldelisle/mutantGenome/"
FourCDirectory="/scratch/ldelisle/4C/"
FourCAnalysis="analysis_TgN3840"

# Check that the sra table is available:
if [ ! -e ${gitHubDirectory}/4C/sraTable.txt ]; then
  echo "sraTable.txt does not exists."
  exit 1
fi

# Check that the mutant genome is available:
if [ ! -e ${mutantGenomeDirectory}/TgN3840 ]; then
  echo "mutant genome does not exists. Please run ${gitHubDirectory}/mutantGenome/prepareMutantGenome.sh first"
  exit 1
fi

# First we need to create all the directories
mkdir -p ${FourCDirectory}

# To prepare the mapseq step:
# Update the init_bein_minilims.py
cp $gitHubDirectory/scripts/4C/init_bein_minilims.py ${FourCDirectory}
sed -i "s#/mnt/BBCF-raw#${FourCDirectory}#g" ${FourCDirectory}/init_bein_minilims.py
sed -i "s#/home/leleu/htsstation#${pathForBBCFutils}#g" ${FourCDirectory}/init_bein_minilims.py

# Begin to launch jobs:

# Get fastq from sra
jidGetFastq=$(sbatch --chdir $FourCDirectory ${gitHubDirectory}/4C/slurm_scripts/01_getDemultiplexedFastq.sh $FourCAnalysis $gitHubDirectory | awk '{print $NF}')

# Launch the mapping
jidMapping4C=$(sbatch --chdir $FourCDirectory --dependency afterok:${jidGetFastq} ${gitHubDirectory}/4C/slurm_scripts/02_mapping.sh $FourCAnalysis $gitHubDirectory $mutantGenomeDirectory | awk '{print $NF}')

# Prepare all files for the 4C
jidpre4C=$(sbatch --chdir $FourCDirectory --dependency afterok:${jidMapping4C} ${gitHubDirectory}/4C/slurm_scripts/03_prepare4C.sh $FourCAnalysis $gitHubDirectory $mutantGenomeDirectory | awk '{print $NF}')

# Launch the 4C
sbatch --chdir $FourCDirectory --dependency afterok:${jidpre4C} ${gitHubDirectory}/4C/slurm_scripts/04_4Cseq.sh $FourCAnalysis $gitHubDirectory
