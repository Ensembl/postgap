#!/bin/bash

mkdir -p databases
cd databases

# GWAS databases
wget https://storage.googleapis.com/postgap-data/GWAS_DB.txt
wget https://storage.googleapis.com/postgap-data/GRASP.txt
wget https://storage.googleapis.com/postgap-data/Phewas_Catalog.txt

# Cisregulatory databases
wget https://storage.googleapis.com/postgap-data/Fantom5.bed.gz
wget https://storage.googleapis.com/postgap-data/Fantom5.bed.gz.tbi
wget https://storage.googleapis.com/postgap-data/Fantom5.fdrs
wget https://storage.googleapis.com/postgap-data/DHS.bed.gz
wget https://storage.googleapis.com/postgap-data/DHS.bed.gz.tbi
wget https://storage.googleapis.com/postgap-data/DHS.fdrs
wget https://storage.googleapis.com/postgap-data/Ensembl_TSSs.bed.gz
wget https://storage.googleapis.com/postgap-data/Ensembl_TSSs.bed.gz.tbi
wget https://storage.googleapis.com/postgap-data/Regulome.bed.gz
wget https://storage.googleapis.com/postgap-data/Regulome.bed.gz.tbi
wget https://storage.googleapis.com/postgap-data/pchic.bed.gz
wget https://storage.googleapis.com/postgap-data/pchic.bed.gz.tbi

# Population structure data
mkdir 1000Genomes
cd 1000Genomes

mkdir CEPH
cd CEPH

wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/CEPH/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi

cd ..

cd ..

cd ..
