#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 10G
#SBATCH --cpus-per-task 10
#SBATCH --time 4:00:00
#SBATCH --array=1-361
#SBATCH --job-name basecall
#SBATCH --chdir /scratch/ldelisle/nCATS/

path="$PWD/"

# We assume all fast5 are in $path/fast5

fast5File=$(ls $path/fast5/*.fast5 | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print}')
uuid=$(uuidgen)
mkdir -p callOnGoing_${uuid}
mkdir -p guppy_output_${uuid}
ln -s $fast5File callOnGoing_${uuid}/
~/softwares/ont-guppy-cpu/bin/guppy_basecaller --version
# : Guppy Basecalling Software, (C) Oxford Nanopore Technologies, Limited. Version 3.1.5+781ed57
~/softwares/ont-guppy-cpu/bin/guppy_basecaller  -i  callOnGoing_${uuid}/ -s guppy_output_${uuid}/ --flowcell FLO-MIN106 --kit SQK-LSK109 --num_callers 1 --cpu_threads_per_caller 28 -r 

mkdir -p guppy_output
cat guppy_output_${uuid}/*.fastq > guppy_output/fastq_${SLURM_ARRAY_TASK_ID}.fastq
rm -r callOnGoing_${uuid}
rm -r guppy_output_${uuid}
