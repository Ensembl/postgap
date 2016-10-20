DEST_DIR=~/lustre2/CTTV24/databases_test

default: download process
download: create_dir d_GRASP d_Phewas_Catalog d_GWAS_DB d_Fantom5 d_DHS d_Regulome d_1000Genomes
process: GRASP Phewas_Catalog GWAS_DB Fantom5 DHS Regulome 1000Genomes

clean_raw:
	rm -rf ${DEST_DIR}/raw/*

clean_all:
	rm -rf ${DEST_DIR}/*

create_dir:
	mkdir -p ${DEST_DIR}
	mkdir -p ${DEST_DIR}/raw

d_GRASP:
	wget -nc https://s3.amazonaws.com/NHLBI_Public/GRASP/GraspFullDataset2.zip -qO ${DEST_DIR}/raw/GRASP.zip

GRASP:
	unzip -c ${DEST_DIR}/raw/GRASP.zip | python scripts/preprocessing/column.py 12 > ${DEST_DIR}/GRASP.txt


d_Phewas_Catalog:
	wget -nc http://phewas.mc.vanderbilt.edu/phewas-catalog.csv > ${DEST_DIR}/raw/Phewas_Catalog.csv

Phewas_Catalog: 
	python scripts/preprocessing/csvToTsv.py ${DEST_DIR}/raw/Phewas_Catalog.csv  > ${DEST_DIR}/Phewas_Catalog.txt

d_GWAS_DB:
	wget -nc ftp://jjwanglab.org/GWASdb/old_release/GWASdb_snp_v4.zip -O ${DEST_DIR}/raw/GWAS_DB.zip
GWAS_DB:
	unzip -c ${DEST_DIR}/raw/GWAS_DB.zip > ${DEST_DIR}/GWAS_DB.txt

GWAS_Catalog:
	wget https://www.ebi.ac.uk/gwas/api/search/downloads/alternative -O ${DEST_DIR}/raw/GWAS_Catalog.txt

d_Fantom5:
	wget -nc http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed -O ${DEST_DIR}/raw/Fantom5.txt

Fantom5:
	cat ${DEST_DIR}/raw/Fantom5.txt | cut -f4 | tr ';' '\t' | cut -f1,3,5 | grep 'FDR:' | sed -e 's/FDR://' -e 's/^chr//' -e 's/-/\t/' -e 's/:/\t/' > ${DEST_DIR}/Fantom5.bed
	cat ${DEST_DIR}/Fantom5.bed | python scripts/preprocessing/STOPGAP_FDR.py > ${DEST_DIR}/Fantom5.fdrs

d_DHS:
	wget -nc ftp://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/openchrom/jan2011/dhs_gene_connectivity/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8.gz -qO ${DEST_DIR}/raw/DHS.txt.gz 

DHS:
	gzip -dc ${DEST_DIR}/raw/DHS.txt.gz | awk 'BEGIN {OFS="\t"} {print $$5,$$6,$$7,$$4,$$8}' | sed -e 's/^chr//' > ${DEST_DIR}/DHS.txt
	cat ${DEST_DIR}/DHS.txt | python scripts/preprocessing/STOPGAP_FDR.py > ${DEST_DIR}/DHS.fdrs

d_Regulome:
	wget -nc http://regulomedb.org/downloads/RegulomeDB.dbSNP132.Category1.txt.gz -qO ${DEST_DIR}/raw/regulome1.csv.gz
	wget -nc http://regulomedb.org/downloads/RegulomeDB.dbSNP132.Category2.txt.gz -qO ${DEST_DIR}/raw/regulome2.csv.gz
	wget -nc http://regulomedb.org/downloads/RegulomeDB.dbSNP132.Category3.txt.gz -qO ${DEST_DIR}/raw/regulome3.csv.gz

Regulome:
	gzip -dc ${DEST_DIR}/raw/regulome1.csv.gz > ${DEST_DIR}/regulome1.csv
	gzip -dc ${DEST_DIR}/raw/regulome2.csv.gz > ${DEST_DIR}/regulome2.csv
	gzip -dc ${DEST_DIR}/raw/regulome3.csv.gz > ${DEST_DIR}/regulome3.csv
	cat ${DEST_DIR}/regulome1.csv ${DEST_DIR}/regulome2.csv ${DEST_DIR}/regulome3.csv > ${DEST_DIR}/regulome.csv
	rm ${DEST_DIR}/regulome1.csv ${DEST_DIR}/regulome2.csv ${DEST_DIR}/regulome3.csv
	awk 'BEGIN {FS="\t"} { print $$1,$$2,$$2 + 1,$$4 }' ${DEST_DIR}/regulome.csv | sed -e 's/^chr//' > ${DEST_DIR}/Regulome.bed
	python scripts/preprocessing/regulome_tidy.py ${DEST_DIR}

d_1000Genomes:
	mkdir -p ./databases/raw/1000Genomes
	cat ./scripts/preprocessing/links.txt | xargs -n1 wget -nc -P ${DEST_DIR}/raw/1000Genomes/

define process_1000Genomes_file
gzip -dc $(1) \
| vcfkeepsamples - `cat ./scripts/preprocessing/CEPH_samples.txt`\
| bcftools convert -Ob \
> ${DEST_DIR}/1000Genomes/CEPH/`basename $(1) | sed -e 's/vcf.gz/bcf.gz/'`;
endef

.PHONY: 1000Genomes

1000Genomes:
	$(eval vcf_files := $(wildcard ${DEST_DIR}/raw/1000Genomes/*.vcf.gz))
	$(foreach file, $(vcf_files), $(call process_1000Genomes_file, $(file)))
	$(eval bcf_files := $(wildcard ${DEST_DIR}/1000Genomes/*.bcf.gz))
	$(foreach file, $(bcf_files), bcftools index $file;)
