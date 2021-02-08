#!/bin/bash -l

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 50G
#SBATCH --cpus-per-task 28
#SBATCH --array=1-50
#SBATCH --time 12:00:00
#SBATCH --job-name 4C

desc=$1
gitHubDirectory=$2

pathWithInstall="$PWD/"
wd="${pathWithInstall}/${desc}"
pathForScripts="${gitHubDirectory}/scripts/4C/"

cd ${wd}

conffile=$(ls $PWD/config_4cseq_*.txt | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $1}')

# We put an array which was on purpose too large
# to be sure to get all config files
# For SLURM_ARRAY_TASK_ID above the number of config files
# We exit with 0
if [ -z ${conffile} ]; then
  echo "There is no conffile"
  if [ ${SLURM_ARRAY_TASK_ID} = "1" ]; then
    echo "There is no conffile at all."
    exit 1
  fi
  exit 0
fi

# The name should be something like invTDOM_split1
name=$(basename ${conffile} .txt | awk '{gsub("config_4cseq_", "", $1); print $1}')


module purge
slmodules -r deprecated
module load gcc/6.4.0
module load expat/2.2.2
module load zlib/1.2.11
module load bowtie2/2.3.4.1
module load picard/2.18.14
module load bowtie/1.2
export PYTHONPATH="/home/ldelisle/softwares/bbcfutils/Python/:/home/ldelisle/lib64/python2.7/site-packages/:/home/ldelisle/lib/python2.7/site-packages/:/home/ldelisle/.local/lib/python2.7/site-packages/:$PYTHONPATH"
export PATH="/home/ldelisle/softwares/bbcfutils/Python/:/home/ldelisle/softwares/bbcfutils/bash_scripts/:/home/ldelisle/softwares/exonerate-2.2.0-x86_64/bin/:/ssoft/spack/paien/v2/opt/spack/linux-rhel7-x86_E5v4_Mellanox/intel-18.0.2/r-3.5.0-py7lua2tnuyszaeerbapgnmjtvmr3jeh/bin/:/home/ldelisle/softwares/samtools-0.1.19/:/home/ldelisle/softwares/FastQC/:$PATH"
export LD_LIBRARY_PATH="/home/ldelisle/bin/lib:$LD_LIBRARY_PATH"

desc=$(cat ${conffile} | awk -F "\047" '$1=="description="{print $2}')
if [ ! -e 4cseq_minilims ]; then
  if [ ! -e ${pathWithInstall}/init_bein_minilims.py ]; then
    echo "${pathWithInstall}/init_bein_minilims.py does not exists"
    exit 1
  fi
  python ${pathWithInstall}/init_bein_minilims.py -d ${wd} -m 4cseq
fi
limspath=${wd}  
run_htsstation.py 4cseq -v local -w ${wd} -c ${conffile} --basepath=${limspath} --no-email
# Rename the files in the minilims to the output folder
if [ ! -e ${pathForScripts}/parse_output.sh ]; then
  echo "${pathForScripts}/parse_output.sh does not exists"
  exit 1
fi
mkdir res_files_4Cseq_${name}/
${pathForScripts}/parse_output.sh 4cseq_minilims ${desc} res_files_4Cseq_${name}/

mkdir -p toGEO
for f in res_files_4Cseq_${name}/*_norm_smoothed_11FragsPerWin.bedGraph.gz; do
  sample=$(basename $f _norm_smoothed_11FragsPerWin.bedGraph.gz | awk '{gsub("segToFrag_", "", $1); print $1}')
  cp $f toGEO/${sample}.bedGraph.gz
  cp ${f/_norm_smoothed_11FragsPerWin.bedGraph.gz/.bw} toGEO/
done
