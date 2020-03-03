#!/bin/bash
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 20
#SBATCH --nodes=1
#SBATCH -J STALLONE
#SBATCH -p vccc
#SBATCH --mem=450G
#SBATCH -o %J.err
#SBATCH -e %J.out
#SBATCH --time=48:00:00
#SBATCH --mail-user=andrew.pattison@unimelb.edu.au
#SBATCH --mail-type=ALL

. /data/cephfs/punim0010/extras/Pattison/miniconda2/etc/profile.d/conda.sh
conda activate STALLONE

snakemake --dag output/final/Somatic_varaints_all_samples.maf output/final/TMB_all_samples.csv | dot -Tsvg > dag.svg
snakemake -j 20 output/final/Somatic_varaints_all_samples.maf output/final/TMB_all_samples.csv
