DEST_DIR=~/hps/postgap/databases
DIR_REGEX=~\/hps\/postgap\/databases

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
	unzip -qc ${DEST_DIR}/raw/GRASP.zip | python preprocessing/pad_columns.py 70 | awk 'BEGIN {FS="\t"} $$11 < 1e-4' | python preprocessing/EFO_suggest.py 13 preprocessing/grasp_suggestions.txt preprocessing/mesh_suggestions.txt preprocessing/GWAS_Catalog_suggestions.txt > ${DEST_DIR}/GRASP.txt

d_Phewas_Catalog:
	wget -nc http://phewascatalog.org/files/phewas-catalog.csv.zip -qO ${DEST_DIR}/raw/Phewas_Catalog.csv.zip
	unzip -d ${DEST_DIR}/raw/ ${DEST_DIR}/raw/Phewas_Catalog.csv.zip

Phewas_Catalog: 
	cat ${DEST_DIR}/raw/phewas-catalog.csv | python preprocessing/csvToTsv.py | python preprocessing/EFO_suggest.py 3 preprocessing/grasp_suggestions.txt preprocessing/mesh_suggestions.txt preprocessing/GWAS_Catalog_suggestions.txt > ${DEST_DIR}/Phewas_Catalog.txt

d_GWAS_DB:
	wget -nc ftp://jjwanglab.org/GWASdb/old_release/GWASdb_snp_v4.zip -qO ${DEST_DIR}/raw/GWAS_DB.zip

GWAS_DB:
	unzip -qc ${DEST_DIR}/raw/GWAS_DB.zip | awk 'BEGIN {FS="\t"; OFS="\t"} NF > 5 {print $$1, $$2, $$3, $$4, $$5, $$6}' | python preprocessing/EFO_suggest.py 6 preprocessing/grasp_suggestions.txt preprocessing/mesh_suggestions.txt preprocessing/GWAS_Catalog_suggestions.txt > ${DEST_DIR}/GWAS_DB.txt

GWAS_Catalog:
	wget https://www.ebi.ac.uk/gwas/api/search/downloads/alternative -qO ${DEST_DIR}/raw/GWAS_Catalog.txt

d_Neale:
	wget https://storage.googleapis.com/postgap-data/postgap_input.nealeUKB_20170915.clumped.1Mb.tsv -q0 ${DEST_DIR}/raw/Neale_UKB.txt

Neale:
	cp ${DEST_DIR}/raw/Neale_UKB.txt ${DEST_DIR}/Neale_UKB.txt

d_Fantom5:
	wget -nc http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed -qO ${DEST_DIR}/raw/Fantom5.txt

Fantom5:
	cat ${DEST_DIR}/raw/Fantom5.txt | cut -f4 | tr ';' '\t' | cut -f1,3,5 | grep 'FDR:' | sed -e 's/FDR://' -e 's/^chr//' -e 's/-/\t/' -e 's/:/\t/' | sort -k1,1 -k2,2n > ${DEST_DIR}/Fantom5.bed
	cat ${DEST_DIR}/Fantom5.bed | python preprocessing/STOPGAP_FDR.py > ${DEST_DIR}/Fantom5.fdrs

d_DHS:
	wget -nc ftp://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/openchrom/jan2011/dhs_transcript_connectivity/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8.gz -qO ${DEST_DIR}/raw/DHS.txt.gz 

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
	wget -q -O - ftp://ftp.ensembl.org/pub/grch37/update/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz | gzip -dc | grep protein_coding | awk 'BEGIN {OFS="\t"} $$3 == "transcript" && $$7== "+" { print $$1, $$4, $$4+1, $$10 } $$3 == "transcript" && $$7== "-" { print $$1, $$5, $$5+1, $$10 } ' | tr -d '";' | sort -k1,1 -k2,2n > ${DEST_DIR}/Ensembl_TSSs.bed

tabix: bgz
	$(eval bgz_files := $(wildcard ${DEST_DIR}/*.bed.gz))
	$(foreach file, $(bgz_files), tabix -f -p bed $(file);)

bgz:
	$(eval bed_files := $(wildcard ${DEST_DIR}/*.bed))
	$(foreach file, $(bed_files), bgzip $(file);)

d_1000Genomes:
	mkdir -p ${DEST_DIR}/raw/1000Genomes
	cat ./preprocessing/links.txt | xargs -n1 wget -nc -P ${DEST_DIR}/raw/1000Genomes/

efo_list:
	cut -f71 ${DEST_DIR}/GRASP.txt | tr ',' '\n' | sort | uniq > ${DEST_DIR}/raw/GRASP.efos.txt
	cut -f10 ${DEST_DIR}/Phewas_Catalog.txt | tr ',' '\n' | sort | uniq > ${DEST_DIR}/raw/Phewas_Catalog.efos.txt
	cut -f7 ${DEST_DIR}/GWAS_DB.txt | tr ',' '\n' | sort | uniq > ${DEST_DIR}/raw/GWAS_DB.efos.txt
	rm -f ${DEST_DIR}/raw/GWAS_Catalog.txt 
	wget -nc http://www.ebi.ac.uk/gwas/api/search/downloads/alternative -qO ${DEST_DIR}/raw/GWAS_Catalog.txt 
	cut -f36 ${DEST_DIR}/raw/GWAS_Catalog.txt | tr -d ' '| tr ',' '\n' | sort | uniq | grep -v MAPPED_TRAIT_URI > ${DEST_DIR}/raw/GWAS_Catalog.efos.txt
	sort -m ${DEST_DIR}/raw/GRASP.efos.txt ${DEST_DIR}/raw/Phewas_Catalog.efos.txt ${DEST_DIR}/raw/GWAS_DB.efos.txt ${DEST_DIR}/raw/GWAS_Catalog.efos.txt | uniq | grep '.' |  grep -v 'N/A' > ${DEST_DIR}/raw/all.efos.txt

script: efo_list
	cat ${DEST_DIR}/raw/all.efos.txt | sed -e 's/.*\///' | sed -e 's/\(.*\)/python postgap_and_tests.py --database_dir ${DIR_REGEX} --efos \1 --output \1.txt/' > all_tests.sh


define process_1000Genomes_file
gzip -dc $(1) \
| vcfkeepsamples - `cat ./preprocessing/EUR_samples.txt`\
| vcftools --vcf - --maf 0.01 --min-alleles 2 --max-alleles 2 --recode --stdout \
| bcftools convert -Ob \
> ${DEST_DIR}/1000Genomes/EUR/`basename $(1) | sed -e 's/vcf.gz/bcf/'`;
endef

.PHONY: 1000Genomes

1000Genomes:
	mkdir -p ${DEST_DIR}/1000Genomes/EUR
	$(eval vcf_files := $(wildcard ${DEST_DIR}/raw/1000Genomes/*.vcf.gz))
	$(foreach file, $(vcf_files), $(call process_1000Genomes_file, $(file)))
	$(eval bcf_files := $(wildcard ${DEST_DIR}/1000Genomes/*.bcf))
	$(foreach file, $(bcf_files), bcftools index $(file);)
