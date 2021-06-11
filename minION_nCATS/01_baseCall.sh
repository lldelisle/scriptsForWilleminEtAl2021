#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 10G
#SBATCH --cpus-per-task 10
#SBATCH --time 48:00:00
#SBATCH --array=0-9
#SBATCH --job-name basecall
#SBATCH --chdir /scratch/ldelisle/nCATS/

path="$PWD/"

# We assume all fast5 are in $path/fast5

fast5Files=$(ls $path/fast5/*${SLURM_ARRAY_TASK_ID}.fast5)
uuid=$(uuidgen)
mkdir -p callOnGoing_${uuid}
mkdir -p guppy_output_${uuid}
ln -s $fast5Files callOnGoing_${uuid}/
~/softwares/ont-guppy-cpu/bin/guppy_basecaller --version
# : Guppy Basecalling Software, (C) Oxford Nanopore Technologies, Limited. Version 3.1.5+781ed57
~/softwares/ont-guppy-cpu/bin/guppy_basecaller  -i  callOnGoing_${uuid}/ \
  -s guppy_output_${uuid}/ --flowcell FLO-MIN106 --kit SQK-LSK109 --num_callers 1 \
  --cpu_threads_per_caller 32 --fast5_out

mkdir -p guppy_output
cat guppy_output_${uuid}/*.fastq > guppy_output/fastq_${SLURM_ARRAY_TASK_ID}.fastq
mkdir -p basecalled_fast5
cp guppy_output_${uuid}/workspace/*fast5 basecalled_fast5/
