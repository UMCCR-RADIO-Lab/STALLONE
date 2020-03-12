#!/bin/bash
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --nodes=1
#SBATCH -J STALLONE
#SBATCH -p vccc
#SBATCH --mem=5G
#SBATCH -o %J.err
#SBATCH -e %J.out
#SBATCH --time=48:00:00
#SBATCH --mail-user=andrew.pattison@unimelb.edu.au
#SBATCH --mail-type=ALL

. /data/cephfs/punim0010/extras/Pattison/miniconda2/etc/profile.d/conda.sh
conda activate STALLONE

# Make a config yaml
rm config.yaml
rm -r logs
mkdir -p logs
echo "samples:" >> config.yaml
inputs=$(ls data/samples/*.R1_001.fastq.gz)

for sample in $inputs
do
        s2=$(basename $sample)
        echo "    ${s2/.R1_001.fastq.gz/}: $sample" >> config.yaml
done

# Build a DAG
snakemake --dag | dot -Tsvg > dag.svg
# Do a dry run
#snakemake -n -j 32 
# Run for real 
snakemake -j 20 --rerun-incomplete --cluster-config cluster.json --cluster "sbatch --time {cluster.time} --cpus-per-task {cluster.cpus-per-task}  -p {cluster.p} --mem {cluster.mem} -o {cluster.o} -e {cluster.e}"

# Cancel all your jobs on spartan (just in case)
#squeue -u andrew_pattison | grep 149 | awk '{print $1}' | xargs -n 1 scancel
