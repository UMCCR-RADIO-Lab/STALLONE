. /data/cephfs/punim0010/extras/Pattison/miniconda2/etc/profile.d/conda.sh
conda activate GDC_downloader

# Need to use the web version of the tool inside the conda env to make it work.
/data/cephfs/punim0648/Pattison_projects/STALLONE/TCGA_manifests/gdc-client download -m /data/cephfs/punim0648/Pattison_projects/STALLONE/TCGA_manifests/SKCM_RNA_manifest.tsv -t /data/cephfs/punim0648/Pattison_projects/STALLONE/TCGA_manifests/gdc-user-token.2020-02-26T02_57_47.831Z.txt \
-d /data/cephfs/punim0648/Pattison_projects/STALLONE/data/TCGA_bams/ -n 10

