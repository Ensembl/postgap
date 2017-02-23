DEST_DIR=databases_dir

default: download process
download: create_dir d_GRASP d_Phewas_Catalog d_GWAS_DB d_Fantom5 d_DHS d_Regulome d_pchic d_1000Genomes
process: GRASP Phewas_Catalog GWAS_DB Fantom5 DHS Regulome tabix pchic 1000Genomes

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
	unzip -qc ${DEST_DIR}/raw/GRASP.zip | python preprocessing/pad_columns.py 70 | awk 'BEGIN {FS="\t"} $$11 < 1e-4' | python preprocessing/EFO_suggest.py 12 > ${DEST_DIR}/GRASP.txt

d_Phewas_Catalog:
	wget -nc http://phewas.mc.vanderbilt.edu/phewas-catalog.csv -qO ${DEST_DIR}/raw/Phewas_Catalog.csv

Phewas_Catalog: 
	cat ${DEST_DIR}/raw/Phewas_Catalog.csv | python preprocessing/csvToTsv.py | python preprocessing/EFO_suggest.py 3 > ${DEST_DIR}/Phewas_Catalog.txt

d_GWAS_DB:
	wget -nc ftp://jjwanglab.org/GWASdb/old_release/GWASdb_snp_v4.zip -qO ${DEST_DIR}/raw/GWAS_DB.zip

GWAS_DB:
	unzip -qc ${DEST_DIR}/raw/GWAS_DB.zip | awk 'BEGIN {FS="\t"; OFS="\t"} NF > 5 {print $$1, $$2, $$3, $$4, $$5, $$6}' | python preprocessing/EFO_suggest.py 6 > ${DEST_DIR}/GWAS_DB.txt

GWAS_Catalog:
	wget https://www.ebi.ac.uk/gwas/api/search/downloads/alternative -qO ${DEST_DIR}/raw/GWAS_Catalog.txt

d_Fantom5:
	wget -nc http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed -qO ${DEST_DIR}/raw/Fantom5.txt

Fantom5:
	cat ${DEST_DIR}/raw/Fantom5.txt | cut -f4 | tr ';' '\t' | cut -f1,3,5 | grep 'FDR:' | sed -e 's/FDR://' -e 's/^chr//' -e 's/-/\t/' -e 's/:/\t/' | sort -k1,1 -k2,2n > ${DEST_DIR}/Fantom5.bed
	cat ${DEST_DIR}/Fantom5.bed | python preprocessing/STOPGAP_FDR.py > ${DEST_DIR}/Fantom5.fdrs

d_DHS:
	wget -nc ftp://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/openchrom/jan2011/dhs_gene_connectivity/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8.gz -qO ${DEST_DIR}/raw/DHS.txt.gz 

DHS:
	gzip -dc ${DEST_DIR}/raw/DHS.txt.gz | awk 'BEGIN {OFS="\t"} {print $$5,$$6,$$7,$$4,$$8}' | sed -e 's/^chr//' | sort -k1,1 -k2,2n > ${DEST_DIR}/DHS.bed
	cat ${DEST_DIR}/DHS.bed | python preprocessing/STOPGAP_FDR.py > ${DEST_DIR}/DHS.fdrs

d_Regulome:
	wget -nc http://regulomedb.org/downloads/RegulomeDB.dbSNP132.Category1.txt.gz -qO ${DEST_DIR}/raw/regulome1.csv.gz
	wget -nc http://regulomedb.org/downloads/RegulomeDB.dbSNP132.Category2.txt.gz -qO ${DEST_DIR}/raw/regulome2.csv.gz
	wget -nc http://regulomedb.org/downloads/RegulomeDB.dbSNP132.Category3.txt.gz -qO ${DEST_DIR}/raw/regulome3.csv.gz

Regulome:
	gzip -dc ${DEST_DIR}/raw/regulome[123].csv.gz | sed -e 's/^chr//' | awk 'BEGIN {FS="\t"; OFS="\t"} { print $$1,$$2,$$2 + 1,$$5 }' | sort -k1,1 -k2,2n > ${DEST_DIR}/Regulome.bed

d_pchic:
	mkdir -p ${DEST_DIR}/raw/pchic
	wget -q -O - ftp://ftp.ebi.ac.uk/pub/contrib/pchic/CHiCAGO/ | tr '<' '\n' | grep '.gz"' | sed -e 's/A HREF="//' -e 's/".*//' | sort | uniq | xargs -n1 -I file wget ftp://ftp.ebi.ac.uk/pub/contrib/pchic/CHiCAGO/file -P ${DEST_DIR}/raw/pchic

pchic: Ensembl
	gzip -dc ${DEST_DIR}/raw/pchic/* | sed -e 's/\<chr//g' | tr ':\-,' '\t' | bedtools intersect -wa -wb -a stdin -b ${DEST_DIR}/Ensembl_TSSs.bed | cut -f4,5,6,13,7 | sort -k1,1 -k2,2n > ${DEST_DIR}/pchic.bed

Ensembl:
	wget -q -O - ftp://ftp.ensembl.org/pub/grch37/update/gtf/homo_sapiens/Homo_sapiens.GRCh37.82.gtf.gz | gzip -dc | grep protein_coding | awk 'BEGIN {OFS="\t"} $$3 == "gene" && $$7== "+" { print $$1, $$4, $$4+1, $$10 } $$3 == "gene" && $$7== "-" { print $$1, $$5, $$5+1, $$10 } ' | tr -d '";' | sort -k1,1 -k2,2n > ${DEST_DIR}/Ensembl_TSSs.bed

tabix: bgz
	$(eval bgz_files := $(wildcard ${DEST_DIR}/*.bed.gz))
	$(foreach file, $(bgz_files), tabix -f -p bed $(file);)

bgz:
	$(eval bed_files := $(wildcard ${DEST_DIR}/*.bed))
	$(foreach file, $(bed_files), bgzip $(file);)

d_1000Genomes:
	mkdir -p ${DEST_DIR}/raw/1000Genomes
	cat ./preprocessing/links.txt | xargs -n1 wget -nc -P ${DEST_DIR}/raw/1000Genomes/

define process_1000Genomes_file
gzip -dc $(1) \
| vcfkeepsamples - `cat ./preprocessing/EUR_samples.txt`\
| bcftools convert -Ob \
> ${DEST_DIR}/1000Genomes/EUR/`basename $(1) | sed -e 's/vcf.gz/bcf/'`;
endef

.PHONY: 1000Genomes

1000Genomes:
	$(eval vcf_files := $(wildcard ${DEST_DIR}/raw/1000Genomes/*.vcf.gz))
	$(foreach file, $(vcf_files), $(call process_1000Genomes_file, $(file)))
	$(eval bcf_files := $(wildcard ${DEST_DIR}/1000Genomes/*.bcf))
	$(foreach file, $(bcf_files), bcftools index $(file);)
