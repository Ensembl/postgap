#!/bin/bash

mkdir -p databases
cd databases

# GWAS databases
wget https://storage.googleapis.com/postgap-data/GWAS_DB.txt
wget https://storage.googleapis.com/postgap-data/GRASP.txt
wget https://storage.googleapis.com/postgap-data/Phewas_Catalog.txt
wget https://storage.googleapis.com/postgap-data/postgap_input.nealeUKB_20170915.clumped.1Mb.tsv -O Neale_UKB.txt


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
wget https://storage.googleapis.com/postgap-data/GERP.bw

# Population structure data
mkdir 1000Genomes
cd 1000Genomes

mkdir AFR
cd AFR

wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AFR/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi

cd ..

mkdir AMR
cd AMR

wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/AMR/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi

cd ..

mkdir EAS
cd EAS

wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EAS/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi

cd ..

mkdir EUR
cd EUR

wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/EUR/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi

cd ..

mkdir SAS
cd SAS

wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf
wget https://storage.googleapis.com/postgap-data/1000Genomes/SAS/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf.csi

cd ..

cd ..

cd ..
