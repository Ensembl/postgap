def GWASCatalog(diseases, efos):
	"""

		Returns all GWAS SNPs associated to a disease in GWAS Catalog
		Args:
		* [ string ] (trait descriptions)
		* [ string ] (trait EFO identifiers)
		Returntype: [ GWAS_Association ]

	"""
	file = open(DATABASES_DIR+"/GWAS_Catalog.txt")
	res = concatenate(get_gwas_catalog_association(line, diseases, efos) for line in file)

	if DEBUG:
		print "\tFound %i GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in GWAS Catalog" % (len(res), ", ".join(diseases), ", ".join(efos))

	return res

def get_gwas_catalog_association(line, diseases, efos):
	'''

		GWAS Catalog flat file format (waiting for REST API...)

		1.	DATE ADDED TO CATALOG: Date added to catalog
		2.	PUBMEDID: PubMed identification number
		3.	FIRST AUTHOR: Last name of first author
		4.	DATE: Publication date (online (epub) date if available)
		5.	JOURNAL: Abbreviated journal name
		6.	LINK: PubMed URL
		7.	STUDY: Title of paper (linked to PubMed abstract)
		8.	DISEASE/TRAIT: Disease or trait examined in study
		9.	INITIAL SAMPLE DESCRIPTION: Sample size for Stage 1 of GWAS
		10.	REPLICATION SAMPLE DESCRIPTION: Sample size for subsequent replication(s)
		11.	REGION: Cytogenetic region associated with rs number (NCBI)
		12.	CHR_ID: Chromosome number associated with rs number (NCBI)
		13.	CHR_POS: Chromosomal position associated with rs number (dbSNP Build 144, Genome Assembly GRCh38.p2, NCBI)
		14.	REPORTED GENE (S): Gene(s) reported by author
		15.	MAPPED GENE(S): Gene(s) mapped to the strongest SNP (NCBI). If the SNP is located within a gene, that gene is listed. If the SNP is intergenic, the upstream and downstream genes are listed, separated by a hyphen.
		16.	UPSTREAM_GENE_ID: Entrez Gene ID for nearest upstream gene to rs number, if not within gene (NCBI)
		17.	DOWNSTREAM_GENE_ID: Entrez Gene ID for nearest downstream gene to rs number, if not within gene (NCBI)
		18.	SNP_GENE_IDS: Entrez Gene ID, if rs number within gene; multiple genes denotes overlapping transcripts (NCBI)
		19.	UPSTREAM_GENE_DISTANCE: distance in kb for nearest upstream gene to rs number, if not within gene (NCBI)
		20.	DOWNSTREAM_GENE_DISTANCE: distance in kb for nearest downstream gene to rs number, if not within gene (NCBI)
		21.	STRONGEST SNP-RISK ALLELE: SNP(s) most strongly associated with trait + risk allele (? for unknown risk allele). May also refer to a haplotype.
		22.	SNPS: Strongest SNP; if a haplotype is reported above, may include more than one rs number (multiple SNPs comprising the haplotype)
		23.	MERGED: denotes whether the SNP has been merged into a subsequent rs record (0 = no; 1 = yes; NCBI) SNP_ID_CURRENT: current rs number (will differ from strongest SNP when merged = 1)
		24.	[Inserted] SNP_ID_CURRENT
		25.	CONTEXT: SNP functional class (NCBI)
		26.	INTERGENIC: denotes whether SNP is in intergenic region (0 = no; 1 = yes; NCBI)
		27.	RISK ALLELE FREQUENCY: Reported risk allele frequency associated with strongest SNP
		28.	P-VALUE: Reported p-value for strongest SNP risk allele (linked to dbGaP Association Browser)
		29.	PVALUE_MLOG: -log(p-value)
		30.	P-VALUE (TEXT): Information describing context of p-value (e.g. females, smokers). Note that p-values are rounded to 1 significant digit (for example, a published pvalue of 4.8 x 10-7 is rounded to 5 x 10-7).
		31.	OR or BETA: Reported odds ratio or beta-coefficient associated with strongest SNP risk allele
		32.	95% CI (TEXT): Reported 95% confidence interval associated with strongest SNP risk allele
		33.	PLATFORM (SNPS PASSING QC): Genotyping platform manufacturer used in Stage 1; also includes notation of pooled DNA study design or imputation of SNPs, where applicable
		34.	[Inserted] CNV
		35.	[Deleted] MAPPED_TRAIT: Mapped Experimental Factor Ontology trait for this study
		36.	MAPPED_TRAIT_URI: URI of the EFO trait

	'''
	items = line.rstrip().split('\t')
	if len(items) != 36:
		print line
		for elem in enumerate(items):
			print "\t".join(map(str, elem))
		assert False, line

	my_efos = re.sub("http://www.ebi.ac.uk/efo/", "", items[35]).split(", ")
	if items[7] in diseases or any(my_efo in efos for my_efo in my_efos):
		return [ 
			GWAS_Association(
				pvalue = float(items[27]),
				snp = snp,
				disease = items[7],
				efo = items[35],
				source = 'GWAS Catalog',
				study = None
			)
			for snp in items[21].split(',')
		]
	else:
		return []

