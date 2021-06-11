#!/bin/bash -l

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 30G
#SBATCH --cpus-per-task 28
#SBATCH --array=1
#SBATCH --time 24:00:00
#SBATCH --job-name mapping

FourCAnalysis=$1
gitHubDirectory=$2
pathWithMutantGenome=$3

pathWithInstall="$PWD/"
# The working directory is where you have all demultiplexed fastqs
wd="${pathWithInstall}/${FourCAnalysis}/"
pathWithTableWithGenomes="${gitHubDirectory}/mutantGenome/table.txt"
pathForScripts="${gitHubDirectory}/scripts/4C/"

# Load all dependencies
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


# Get the genome and brFiles from the table
if [ ! -e ${pathWithTableWithGenomes} ]; then
  echo "${pathWithTableWithGenomes} does not exists"
  exit 1
fi

genome=$(cat ${pathWithTableWithGenomes} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $1}')
brFiles=$(cat ${pathWithTableWithGenomes} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $2}')

# Attribute the good name
if [ -z ${brFiles} ]; then
  assembly=${genome}
  if [ ${genome} = "mm10" ]; then
    genome="Wt"
  fi
else
  assembly=mm10_${genome}
fi

mkdir -p ${pathWithInstall}/genrep/nr_assemblies/fasta/
mkdir -p ${pathWithInstall}/genrep/nr_assemblies/bowtie2/
mkdir -p ${pathWithInstall}/genrep/nr_assemblies/info_json/

# Go do the working directory
cd ${wd}
# Generate the config file
if [ ! -e ${pathForScripts}/build_config_mapseq.sh ]; then
  echo "${pathForScripts}/build_config_mapseq.sh does not exists"
  exit 1
fi
# The build_config_mapseq.sh expect full path
bash ${pathForScripts}/build_config_mapseq.sh $PWD/${genome}/*.fastq > config_mapseq_${SLURM_ARRAY_TASK_ID}.txt

# Adapt it to the genome
sed -i "s/'mappingAfter/'mapping_${SLURM_ARRAY_TASK_ID}_/g" config_mapseq_${SLURM_ARRAY_TASK_ID}.txt
sed -i "s/assembly_id='mm10'/assembly_id='${assembly}'/g" config_mapseq_${SLURM_ARRAY_TASK_ID}.txt

# If it is not in the genrep directory then add it
if [ ! -e ${pathWithInstall}/genrep/nr_assemblies/fasta/${assembly}.fa.gz ]; then
  if [ ! -e ${pathWithMutantGenome}/${genome}/${assembly}.fa.gz ]; then
    if [ ! -e ${pathWithMutantGenome}/${genome}/${assembly}.fa ]; then
      echo "There is no ${pathWithMutantGenome}/${genome}/${assembly}.fa"
      exit 1
    fi
    gzip -c ${pathWithMutantGenome}/${genome}/${assembly}.fa > ${pathWithMutantGenome}/${genome}/${assembly}.fa.gz
  fi
  cp ${pathWithMutantGenome}/${genome}/${assembly}.fa.gz ${pathWithInstall}/genrep/nr_assemblies/fasta/
  if [ ! -e ${pathWithMutantGenome}/${genome}/${assembly}.1.bt2 ]; then
    echo "${pathWithMutantGenome}/${genome}/${assembly}.1.bt2 does not exists"
    exit 1
  fi
  cp ${pathWithMutantGenome}/${genome}/${assembly}.*.bt2 ${pathWithInstall}/genrep/nr_assemblies/bowtie2/
  if [ ! -e ${pathWithMutantGenome}/${genome}/${assembly}.json ]; then
    echo "${pathWithMutantGenome}/${genome}/${assembly}.json does not exists"
    exit 1
  fi
  cp ${pathWithMutantGenome}/${genome}/${assembly}.json ${pathWithInstall}/genrep/nr_assemblies/info_json/
  if [ ! -e ${pathWithMutantGenome}/${genome}/${assembly}_chromosomes.json ]; then
    echo "${pathWithMutantGenome}/${genome}/${assembly}_chromosomes.json does not exists"
    exit 1
  fi
  cp ${pathWithMutantGenome}/${genome}/${assembly}_chromosomes.json ${pathWithInstall}/genrep/nr_assemblies/info_json/
fi

# Create a directory for the results
mkdir -p res_files_mapping_${genome}

# Create a minilims for mapping
if [ ! -e mapseq_minilims ]; then
  if [ ! -e ${pathWithInstall}/init_bein_minilims.py ]; then
    echo "${pathWithInstall}/init_bein_minilims.py does not exists"
    exit 1
  fi
  python ${pathWithInstall}/init_bein_minilims.py -d ${wd} -m mapseq
fi
limspath=${wd}

# Run the mapseq part of the pipeline
run_htsstation.py mapseq -v local -w ${wd} -c config_mapseq_${SLURM_ARRAY_TASK_ID}.txt --basepath=${limspath} --no-email

# Get the description from the config file
desc=$(cat config_mapseq_${SLURM_ARRAY_TASK_ID}.txt | awk -F "\047" '$1=="description="{print $2}')

# Rename the files in the minilims to the output folder
if [ ! -e ${pathForScripts}/parse_output.sh ]; then
  echo "${pathForScripts}/parse_output.sh does not exists"
  exit 1
fi
${pathForScripts}/parse_output.sh mapseq_minilims ${desc} res_files_mapping_${genome}

# Extract mapping rate:
if [ ! -e ${pathForScripts}/getMappingRate.py ]; then
  echo "${pathForScripts}/getMappingRate.py does not exists"
  exit 1
fi
python ${pathForScripts}/getMappingRate.py --directory res_files_mapping_${genome} \
  --output res_files_mapping_${genome}/mapping_summary.txt
