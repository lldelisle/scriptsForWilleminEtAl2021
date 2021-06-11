#!/bin/bash

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 100G
#SBATCH --cpus-per-task 36
#SBATCH --time 48:00:00
#SBATCH --job-name combineFast5
#SBATCH --chdir /scratch/ldelisle/nCATS/

path="$PWD/"
module purge
module load gcc/7.4.0
module load python/3.7.3

# The ont-fast5-api has been installed through pip
# pip install --user ont-fast5-api
# We assume all fast5 to concatenate are in $path/basecalled_fast5/

mkdir -p basecalled_fast5_single
multi_to_single_fast5 -i basecalled_fast5/ -s basecalled_fast5_single/ -t 36

mkdir -p basecalled_fast5_combined
single_to_multi_fast5 -i basecalled_fast5_single/ -s basecalled_fast5_combined/ -t 36 -f TgN3840_nCATS -n 1000000 --recursive
