#!/bin/bash
#SBATCH --ntasks-per-node 32
#SBATCH --cpus-per-task 1
#SBATCH --nodes=1
#SBATCH -J VCF2MAF
#SBATCH -p vccc
#SBATCH --mem=300G
#SBATCH -o %J.err
#SBATCH -e %J.out
#SBATCH --time=120:00:00
#SBATCH --mail-user=andrew.pattison@unimelb.edu.au
#SBATCH --mail-type=ALL


# Take the early resuls of the GATK haplotype caller and run vcf2maf
# I can then read them into R for some AF based filtering

vcf2maf(){

vcf=$1
output=$2

maf=$(basename $vcf)
maf=$(echo $maf | cut -d '.' -f 1)
maf_out=$output/$maf.maf

ref=/data/cephfs/punim0010/local/development/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa
filt=/data/cephfs/punim0010/local/development.broken/bcbio/genomes/Hsapiens/GRCh37/variation/exac.vcf.gz
vep_path=/data/cephfs/punim0010/extras/Pattison/miniconda/envs/maftools/bin/
vep_data=/data/cephfs/punim0648/Pattison_projects/A5/BCBio_VCFs/VEP_GRCh37_cache/

vcf2maf.pl --input-vcf $vcf --output-maf $maf_out --ref-fasta $ref --filter-vcf $filt \
  --vep-path $vep_path --vep-data $vep_data --vep-forks 8 --tumor-id $maf \
  --species homo_sapiens --retain-info DB,MQ,QD,FS,AD,DP

}

export -f vcf2maf

# Remove any old VEP files from previous runs
rm /data/cephfs/punim0648/Pattison_projects/TCGA-RNA-Seq/maftools/raw_vcfs/*vep*

. /data/cephfs/punim0010/extras/Pattison/miniconda2/etc/profile.d/conda.sh
conda activate vcf2maf

# Prior to running this script
vcfs=$(ls /data/cephfs/punim0648/Pattison_projects/TCGA-RNA-Seq/maftools/raw_vcfs/*.vcf)
outdir=/data/cephfs/punim0648/Pattison_projects/TCGA-RNA-Seq/maftools/raw_mafs

# Run with a whole node

parallel --jobs 4 vcf2maf ::: $vcfs ::: $outdir
