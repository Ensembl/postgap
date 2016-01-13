DEST_DIR=~/lustre2/CTTV24/databases

default: GRASP Phewas_Catalog GWAS_Catalog GWAS_DB Fantom5 DHS Regulome Phenotypes

GRASP: ${DEST_DIR}/GRASP.txt
Phewas_Catalog: ${DEST_DIR}/Phewas_Catalog.txt
GWAS_DB: ${DEST_DIR}/GWAS_DB.txt
GWAS_Catalog: ${DEST_DIR}/GWAS_Catalog.txt
Fantom5: ${DEST_DIR}/Fantom5.txt
DHS: ${DEST_DIR}/DHS.txt
Regulome: ${DEST_DIR}/Regulome.txt
Phenotypes: ${DEST_DIR}/Phenotypes.txt

${DEST_DIR}/GRASP.txt: ${DEST_DIR}
	wget https://s3.amazonaws.com/NHLBI_Public/GRASP/GraspFullDataset2.zip -O ${DEST_DIR}/GRASP.zip
	unzip -c ${DEST_DIR}/GRASP.zip | python scripts/preprocessing.py 12 > ${DEST_DIR}/GRASP.txt

${DEST_DIR}/Phewas_Catalog.txt: ${DEST_DIR}
	wget http://phewas.mc.vanderbilt.edu/phewas-catalog.csv -qO- | python csvToTsv > ${DEST_DIR}/Phewas_Catalog.txt

${DEST_DIR}/GWAS_DB.txt: ${DEST_DIR}
	wget http://jjwanglab.org:8080/gwasdb/GWASdb_ld_snp_v4.zip -O ${DEST_DIR}/GWAS_DB.zip
	unzip -c ${DEST_DIR}/GWAS_DB.zip > ${DEST_DIR}/GWAS_DB.txt

${DEST_DIR}/GWAS_Catalog.txt: ${DEST_DIR}
	wget https://www.ebi.ac.uk/gwas/api/search/downloads/alternative -O ${DEST_DIR}/GWAS_Catalog.txt

${DEST_DIR}/Fantom5.txt: ${DEST_DIR}
	wget http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed -O ${DEST_DIR}/Fantom5.txt
	cat ${DEST_DIR}/Fantom5.txt | perl STOPGAP_FDR.pl | sed -e 's/^chr//' > ${DEST_DIR}/Fantom5.fdrs

${DEST_DIR}/DHS.txt: ${DEST_DIR}
	wget ftp://ftp.ebi.ac.uk/pub/${DEST_DIR}/ensembl/encode/integration_data_jan2011/byDataType/openchrom/jan2011/dhs_gene_connectivity/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8.gz -qO- | gzip -dc | awk 'BEGIN {OFS="\t"} {print $5,$6,$7,$4,$8}' | sed -e 's/^chr//' > ${DEST_DIR}/DHS.txt
	cat ${DEST_DIR}/DHS.txt | perl STOPGAP_FDR.pl > ${DEST_DIR}/DHS.fdrs

${DEST_DIR}/Regulome.txt: ${DEST_DIR}
	wget http://regulomedb.org/downloads/RegulomeDB.dbSNP132.Category1.txt.gz -qO- | gzip -dc > ${DEST_DIR}/regulome.tmp
	wget http://regulomedb.org/downloads/RegulomeDB.dbSNP132.Category2.txt.gz -qO- | gzip -dc >> ${DEST_DIR}/regulome.tmp
	wget http://regulomedb.org/downloads/RegulomeDB.dbSNP132.Category3.txt.gz -qO- | gzip -dc >> ${DEST_DIR}/regulome.tmp
	awk 'BEGIN {FS="\t"} { print $1,$2,$2 + 1,$4 }' ${DEST_DIR}/regulome.tmp | sed -e 's/^chr//' > ${DEST_DIR}/Regulome.txt

${DEST_DIR}/Phenotypes.txt: ${DEST_DIR}:

${DEST_DIR}:
	mkdir -p ${DEST_DIR}
