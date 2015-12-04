default: GRASP Phewas_Catalog GWAS_Catalog GWAS_DB Fantom5

GRASP: databases/GRASP.txt
Phewas_Catalog: databases/Phewas_Catalog.txt
GWAS_DB: databases/GWAS_DB.txt
GWAS_Catalog: databases/GWAS_Catalog.txt
Fantom5: databases/Fantom5.txt

databases/GRASP.txt: databases
	curl https://s3.amazonaws.com/NHLBI_Public/GRASP/GraspFullDataset2.zip > databases/GRASP.zip
	unzip -c databases/GRASP.zip > databases/GRASP.txt

databases/Phewas_Catalog.txt: databases
	curl http://phewas.mc.vanderbilt.edu/phewas-catalog.csv | python csvToTsv > databases/Phewas_Catalog.txt

databases/GWAS_DB.txt: databases
	curl http://jjwanglab.org:8080/gwasdb/GWASdb_ld_snp_v4.zip > databases/GWAS_DB.zip
	unzip -c databases/GWAS_DB.zip > databases/GWAS_DB.txt

databases/GWAS_Catalog.txt: databases
	curl https://www.ebi.ac.uk/gwas/api/search/downloads/full > databases/GWAS_Catalog.txt

databases/Fantom5.txt: databases
	curl http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed > databases/Fantom5.txt

databases:
	mkdir -p databases
  
