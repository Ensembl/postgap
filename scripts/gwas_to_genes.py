#! /usr/bin/env python

"""

Copyright [1999-2016] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License")
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

		 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

"""

"""

	Please email comments or questions to the public Ensembl
	developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

	Questions may also be sent to the Ensembl help desk at
	<http://www.ensembl.org/Help/Contact>.

"""
import os
import sys
import argparse
import re
import requests
import collections
import xmltodict
import pybedtools
import math
import json
import cPickle as pickle
import time
import subprocess

# Class definitions
SNP = collections.namedtuple("SNP", ['rsID', 'chrom', 'pos'])
Gene = collections.namedtuple("Gene", ['name', 'id', 'chrom', 'tss', 'biotype'])
Disease = collections.namedtuple('Disease', ['name', 'efo'])
GWAS_Association = collections.namedtuple('GWAS_Association', ['snp','disease','pvalue','source','study'])
GWAS_SNP = collections.namedtuple('GWAS_SNP', ['snp','pvalue', 'evidence'])
GWAS_Cluster = collections.namedtuple('GWAS_Cluster', ['gwas_snps','ld_snps'])
Cisregulatory_Evidence = collections.namedtuple('Cisregulatory_Evidence', ['snp','gene','score','source','study','tissue'])
Regulatory_Evidence = collections.namedtuple('Regulatory_Evidence', ['snp','score','source','study','tissue'])
GeneSNP_Association = collections.namedtuple('GeneSNP_Association', ['gene', 'snp', 'score', 'cisregulatory_evidence', 'regulatory_evidence'])
GeneCluster_Association = collections.namedtuple('GeneCluster_Association', ['gene', 'cluster', 'score', 'evidence'])
FDR_Model = collections.namedtuple('FDR_Model', ['FDR','BIN_WIDTH','MAX_DISTANCE'])

# Bit and bobs
phenotype_cache = ()
known_genes = {}
VEP_impact_to_score = {
	'HIGH': 4,
	'MEDIUM': 3,
	'LOW': 2,
	'MODIFIER': 1,
	'MODERATE': 1
}

NCBI_Taxon_ID = {
	'Human': 9606
}

# Globals
DATABASES_DIR = None
SPECIES = None
DEBUG = True
PVALUE_CUTOFF = 1e-4

'''
Development TODO list:

A. Datasets to integrate:
	-Cisregulatory annotations:
		PCHIC (STOPGAP Scoring: Single cell line: +1, multiple cell lines: +2)

	-Epigenetic activity:
		PhyloP (STOPGAP Scoring: FPR 0-0.6: +2, 0.6-0.85: +1,0.85-1: +0)

B. Code improvements:
	Pathways analysis (Downstream)
	Take into account population composition in LD calcs

C. Model improvements:
	Replace PICS with Bayesian model
	Fine mapping of summary data
	Tissue selection

'''

def main():
	"""

		Reads commandline parameters, prints corresponding associated genes with evidence info

	"""
	options = get_options()
	if options.output is None:
		output = sys.stdout
	else:
		output = open(options.output, "w")

	if options.rsID is None:
		res = diseases_to_genes(options.diseases, options.efos, "CEPH", options.tissues)
	else:
		res = rsIDs_to_genes(options.rsID, options.tissues)

	if options.json_output:
		formatted_results = "\n".join(map(json.dumps, res))
	elif options.rsID is None:
		formatted_results = pretty_output(res)
	else:
		formatted_results = pretty_snp_output(res)

	output.write(formatted_results)

def get_options():
    """

        Reads commandline parameters
        Returntype:
            {
                diseases: [ string ],
                populations: [ string ],
                tissues: [ string ],
            }

    """
    parser = argparse.ArgumentParser(description='Search GWAS/Regulatory/Cis-regulatory databases for causal genes')
    parser.add_argument('--efos', nargs='*')
    parser.add_argument('--diseases', nargs='*')
    parser.add_argument('--rsID')
    # parser.add_argument('--populations', nargs='*', default=['1000GENOMES:phase_3:GBR'])
    parser.add_argument('--tissues', nargs='*')
    parser.add_argument('--output')
    parser.add_argument('--species', nargs='*', default = 'Human')
    parser.add_argument('--database_dir', dest = 'databases', default = 'databases')
    parser.add_argument('--debug', '-g', action = 'store_true')
    parser.add_argument('--json_output', '-j', action = 'store_true')
    options = parser.parse_args()

    global DATABASES_DIR
    DATABASES_DIR = options.databases
    global SPECIES
    SPECIES = options.species
    global DEBUG
    DEBUG = DEBUG or options.debug

    assert DATABASES_DIR is not None
    assert options.rsID is None or (options.efos is None and options.diseases is None)
    assert options.rsID is not None or options.efos is not None or options.diseases is not None

    if options.diseases is None:
        options.diseases = []

    if options.efos is None:
        options.efos = filter(lambda X: X is not None, (efo_suggest(disease) for disease in options.diseases))

    # Expand list of EFOs to children, concatenate, remove duplicates
    options.efos = concatenate(map(efo_children, options.efos))

    return options

def efo_suggest(term):
	"""

		Find most appropriate EFO term for arbitrary string
		Arg:
		* string
		Returntype: string (EFO ID)

	"""
	server = 'http://www.ebi.ac.uk/spot/zooma/v2/api'
	url_term = re.sub(" ", "%20", term)
	ext = "/summaries?query=%s" % (url_term)
	result = get_rest_json(server, ext)
	'''

		Example search result:
		{
			id: "864A7A7335109CD59C2986398637CB519F21DB05",
			semanticTags: [
			"http://www.orpha.net/ORDO/Orphanet_308410"
			],
			annotationURIs: [
					"http://rdf.ebi.ac.uk/resource/zooma/efo/DE1004D39B114BDEEFB7A2EA93D881AD"
					],
			annotationSourceURIs: [
					"http://www.ebi.ac.uk/efo/efo.owl"
					],
			annotationSummaryTypeName: "disease; Orphanet_308410",
			annotatedPropertyType: "disease",
			annotatedPropertyValue: "Autism-epilepsy syndrome due to branched chain ketoacid dehydrogenase kinase deficiency",
			annotatedPropertyUri: "http://rdf.ebi.ac.uk/resource/zooma/FB9E6464F1830B24C08204157E7A088E",
			quality: 72.82359,
			uri: "http://rdf.ebi.ac.uk/resource/zooma/annotation_summary/864A7A7335109CD59C2986398637CB519F21DB05"
		},


	'''

	hits = filter(lambda X: len(X['semanticTags']) == 1, result)
	if len(hits):
		sorted_hits = sorted(hits, key = lambda X: X['quality'])
		selection = sorted_hits[-1]['semanticTags'][0]
		efo = re.sub('.*/', '\1', selection)
		if DEBUG:
			print "Suggested EFO %s" % efo
		return efo

	else:
		if DEBUG:
			print "Suggested EFO N/A"
		return None

def efo_children(efo):
	"""

		Return list of children EFO IDs
		Arg:
		* string (EFO ID)
		Returntype: [ string ] (EFI IDs)

	"""
	server = 'http://www.ebi.ac.uk'
	page = 1
	res = [efo]

	while (True):
		ext = "/ols/beta/api/ontologies/efo/terms/http253A252F252Fwww.ebi.ac.uk252Fefo252F%s/descendants?page%i&size=10" % (efo, page)
		hash = get_rest_json(server, ext)
		'''

			OWL Output format:
			{
				"_links" : {
					"first" : {
			"href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http253A252F252Fwww.ebi.ac.uk252Fefo252FEFO_0001071/descendants?page=0&size=10"
					},
					"prev" : {
			"href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http253A252F252Fwww.ebi.ac.uk252Fefo252FEFO_0001071/descendants?page=0&size=10"
					},
					"self" : {
			"href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http253A252F252Fwww.ebi.ac.uk252Fefo252FEFO_0001071/descendants?page=1&size=10"
					},
					"last" : {
			"href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http253A252F252Fwww.ebi.ac.uk252Fefo252FEFO_0001071/descendants?page=1&size=10"
					}
				},
				"_embedded" : {
					"terms" : [ {
						"iri" : "http://www.ebi.ac.uk/efo/EFO_1000333",
						"label" : "Lung Inflammatory Myofibroblastic Tumor",
						"description" : [ "An intermediate fibroblastic neoplasm arising from the lung. It is characterized by the presence of spindle-shaped fibroblasts and myofibroblasts, and a chronic inflammatory infiltrate composed of eosinophils, lymphocytes and plasma cells." ],
						"annotation" : {
							"NCI_Thesaurus_definition_citation" : [ "NCIt:C39740" ]
						},
						"synonyms" : null,
						"ontology_name" : "efo",
						"ontology_prefix" : "EFO",
						"ontology_iri" : "http://www.ebi.ac.uk/efo",
						"is_obsolete" : false,
						"is_defining_ontology" : true,
						"has_children" : false,
						"is_root" : false,
						"short_form" : "EFO_1000333",
						"obo_id" : "EFO:1000333",
						"_links" : {
							"self" : {
								"href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http253A252F252Fwww.ebi.ac.uk252Fefo252FEFO_1000333"
							},
							"parents" : {
								"href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http253A252F252Fwww.ebi.ac.uk252Fefo252FEFO_1000333/parents"
							},
							"ancestors" : {
								"href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http253A252F252Fwww.ebi.ac.uk252Fefo252FEFO_1000333/ancestors"
							},
							"jstree" : {
								"href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http253A252F252Fwww.ebi.ac.uk252Fefo252FEFO_1000333/jstree"
							},
							"graph" : {
								"href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http253A252F252Fwww.ebi.ac.uk252Fefo252FEFO_1000333/graph"
							}
						}
					} ]
				},
				"page" : {
					"size" : 10,
					"totalElements" : 20,
					"totalPages" : 2,
					"number" : 1
				}
			}

		'''

		if '_embedded' in hash:
			res += [ result['short_form'] for result in hash['_embedded']['terms'] ]

		if page > int(hash['page']['totalPages']):
			break

		page += 1

	if DEBUG:
		print "EFO children: " + "\t".join(res)
	return res

def diseases_to_genes(diseases, efos, populations, tissues):
	"""

		Associates genes from a list of diseases
		Args:
		* [ string ] (trait descriptions - free strings)
		* [ string ] (trait EFO identifiers)
		* [ string ] (population names)
		* [ string ] (tissue names)
		Returntype: [ GeneCluster_Association ]

	"""
	return gwas_snps_to_genes(diseases_to_gwas_snps(diseases, efos), populations, tissues)

def diseases_to_gwas_snps(diseases, efos):
	"""

		Associates gwas_snps from a list of diseases
		Args:
		* [ string ] (trait descriptions - free strings )
		* [ string ] (trait EFO identifiers)
		Returntype: [ Association_SNP ]

	"""
	res = filter(lambda X: X.pvalue < PVALUE_CUTOFF, scan_disease_databases(diseases, efos))

	if DEBUG:
		print "Found %i GWAS SNPs associated to diseases (%s) or EFO IDs (%s) after p-value filter (%f)" % (len(res), ", ".join(diseases), ", ".join(efos), PVALUE_CUTOFF)

	return res

def scan_disease_databases(diseases, efos):
	"""

		Associates gwas_snps from a list of diseases
		Args:
		* [ string ] (trait descriptions)
		* [ string ] (trait EFO identifiers)
		Returntype: [ GWAS_SNP ]

	"""
	if DEBUG:
		print "Searching for GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in all databases" % (", ".join(diseases), ", ".join(efos))

	hits = concatenate(function(diseases, efos) for function in database_functions)
	associations_by_snp = dict()
	for hit in hits:
		if hit.snp in associations_by_snp:
			record = associations_by_snp[hit.snp]
			record.evidence.append(hit)
			if record.pvalue > hit.pvalue:
				associations_by_snp[hit.snp] = GWAS_SNP(
					snp = record.snp,
					pvalue = hit.pvalue,
					evidence = record.evidence
				)
		else:
			associations_by_snp[hit.snp] = GWAS_SNP(
				snp = hit.snp,
				pvalue = hit.pvalue,
				evidence = [ hit ]
			)

	res = associations_by_snp.values()

	if DEBUG:
		print "Found %i unique GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in all databases" % (len(res), ", ".join(diseases), ", ".join(efos))

	return res

def GWASCatalog(diseases, efos):
	"""

		Returns all GWAS SNPs associated to a disease in GWAS Catalog
		Args:
		* [ string ] (trait descriptions)
		* [ string ] (trait EFO identifiers)
		Returntype: [ GWAS_Association ]

	"""
	if efos is not None and len(efos) > 0:
		res = concatenate(query_gwas_catalog(query) for query in efos)
	else:
		res = concatenate(query_gwas_catalog(query) for query in diseases)

	if DEBUG:
		print "\tFound %i GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in GWAS Catalog" % (len(res), ", ".join(diseases), ", ".join(efos))

	return res

def query_gwas_catalog(term):
	server = 'http://www.ebi.ac.uk'
	url_term = re.sub(" ", "%20", term)
	ext1 = '/gwas/api/search/moreresults?q="%s"&max=0&facet=association&pvalfilter=&orfilter=&betafilter=&datefilter=&sort=' % (url_term)
	hash = get_rest_json(server, ext1)
	try:
		count = int(hash['response']['numFound'])
	except:
		print "Failed on %s%s" % (server, ext1)
	ext2 = '/gwas/api/search/moreresults?q="%s"&max=%i&facet=association&pvalfilter=&orfilter=&betafilter=&datefilter=&sort=' % (url_term, count)
	hash = get_rest_json(server, ext2)
	"""
		{
		  "_version_": 1526208365253886000,
		  "resourcename": "association",
		  "id": "association:21946",
		  "parent": [
		    "disposition",
		  ],
		  "efoLink": [
		    "type II diabetes mellitus|EFO_0001360|http://www.ebi.ac.uk/efo/EFO_0001360"
		  ],
		  "synonym_autosuggest_e": [
		    "T2DM - Type 2 Diabetes mellitus",
		  ],
		  "synonym_autosuggest_ws": [
		    "T2DM - Type 2 Diabetes mellitus",
		  ],
		  "synonym_autosuggest": [
		    "T2DM - Type 2 Diabetes mellitus",
		  ],
		  "synonym": [
		    "T2DM - Type 2 Diabetes mellitus",
		  ],
		  "label_autosuggest_e": [
		    "type II diabetes mellitus"
		  ],
		  "label_autosuggest_ws": [
		    "type II diabetes mellitus"
		  ],
		  "label_autosuggest": [
		    "type II diabetes mellitus"
		  ],
		  "label": [
		    "type II diabetes mellitus"
		  ],
		  "shortform_autosuggest": [
		    "EFO_0001360"
		  ],
		  "shortForm": [
		    "EFO_0001360"
		  ],
		  "traitUri": [
		    "http://www.ebi.ac.uk/efo/EFO_0001360"
		  ],
		  "mappedUri": [
		    "http://www.ebi.ac.uk/efo/EFO_0001360"
		  ],
		  "mappedLabel": [
		    "type ii diabetes mellitus"
		  ],
		  "traitName_s": "Type 2 diabetes",
		  "traitName": [
		    "Type 2 diabetes"
		  ],
		  "ancestryLinks": [
		    "initial|NR|NR|60647",
		  ],
		  "chromosomeName": [
		    "6"
		  ],
		  "studyId": "6787",
		  "reportedGene": [
		    "HLA-DQA2"
		  ],
		  "mappedGeneLinks": [
		    "LOC102725019|102725019"
		  ],
		  "mappedGene": [
		    "LOC102725019"
		  ],
		  "region": [
		    "6p21.32"
		  ],
		  "context": [
		    "nearGene-5"
		  ],
		  "strongestAllele": [
		    "rs3916765-A"
		  ],
		  "riskFrequency": "0.12",
		  "qualifier": [
		    "(Lean)"
		  ],
		  "pValueMantissa": 1,
		  "pValueExponent": -6,
		  "orPerCopyNum": 1.21,
		  "orPerCopyRange": "[1.12-1.31]",
		  "orType": true,
		  "rsId": [
		    "rs3916765"
		  ],
		  "chromosomePosition": [
		    32717773
		  ],
		  "locusDescription": "Single variant",
		  "pubmedId": "22693455",
		  "title": "Stratifying type 2 diabetes cases by BMI identifies genetic risk variants in LAMA1 and enrichment for risk variants in lean compared to obese cases.",
		  "author": [
		    "Perry JR"
		  ],
		  "author_s": "Perry JR",
		  "publication": "PLoS Genet",
		  "publicationDate": "2012-05-30T23:00:00Z",
		  "catalogPublishDate": "2012-08-07T04:10:19Z",
		  "publicationLink": [
		    "Perry JR|2012|22693455"
		  ],
		  "platform": "NR",
		  "initialSampleDescription": "2,112 lean type 2 diabetes cases, 4,123 obese type 2 diabetes cases, 54,412 controls",
		  "replicateSampleDescription": "2,881 lean type 2 diabetes cases, 8,702 obese type 2 diabetes cases, 18,957 controls",
		  "ancestralGroups": [
		    "NR"
		  ],
		  "countriesOfRecruitment": [
		    "NR"
		  ],
		  "numberOfIndividuals": [
		    60647,
		    30540
		  ]
		}


	"""
	try:
		hits = hash['response']['docs']
		return [
			GWAS_Association(
				pvalue = float(hit['pValueMantissa']) * 10**float(hit['pValueExponent']),
				snp = hit['rsId'][0].strip(),
				disease = Disease(name = hit['traitName_s'], efo = hit['shortForm'][0]),
				source = "GWAS Catalog",
				study = "PMID" + hit['pubmedId']
			)
			for hit in hits
		]
	except:
		print "Failed on %s%s" % (server, ext)

def GRASP(diseases, efos):
	"""

		Returns all GWAS SNPs associated to a disease in GRASP
		Args:
		* [ string ] (trait descriptions)
		* [ string ] (trait EFO identifiers)
		Returntype: [ GWAS_Association ]

	"""
	file = open(DATABASES_DIR+"/GRASP.txt")
	res = [ get_grasp_association(line, diseases, efos) for line in file ]
	res = filter(lambda X: X is not None, res)

	if DEBUG:
		print "\tFound %i GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in GRASP" % (len(res), ", ".join(diseases), ", ".join(efos))

	return res

def get_grasp_association(line, diseases, efos):
	'''

		GRASP file format:
		1. NHLBIkey
		2. HUPfield
		3. LastCurationDate
		4. CreationDate
		5. SNPid(dbSNP134)
		6. chr(hg19)
		7. pos(hg19)
		8. PMID
		9. SNPid(in paper)
		10. LocationWithinPaper
		11. Pvalue
		12. Phenotype
		13. PaperPhenotypeDescription
		14. PaperPhenotypeCategories
		15. DatePub
		16. InNHGRIcat(as of 3/31/12)
		17. Journal
		18. Title
		19. IncludesMale/Female Only Analyses
		20. Exclusively Male/Female
		21. Initial Sample Description
		22. Replication Sample Description
		23. Platform [SNPs passing QC]
		24. GWASancestryDescription
		25. TotalSamples(discovery+replication)
		26. TotalDiscoverySamples
		27. European Discovery
		28. African Discovery
		29. East Asian Discovery
		30. Indian/South Asian Discovery
		31. Hispanic Discovery
		32. Native Discovery
		33. Micronesian Discovery
		34. Arab/ME Discovery
		35. Mixed Discovery
		36. Unspecified Discovery
		37. Filipino Discovery
		38. Indonesian Discovery
		39. Total replication samples
		40. European Replication
		41. African Replication
		42. East Asian Replication
		43. Indian/South Asian Replication
		44. Hispanic Replication
		45. Native Replication
		46. Micronesian Replication
		47. Arab/ME Replication
		48. Mixed Replication
		49. Unspecified Replication
		50. Filipino Replication
		51. Indonesian Replication
		52. InGene
		53. NearestGene
		54. InLincRNA
		55. InMiRNA
		56. InMiRNABS
		57. dbSNPfxn
		58. dbSNPMAF
		59. dbSNPalleles/het/se
		60. dbSNPvalidation
		XX61. dbSNPClinStatus
		XX62. ORegAnno
		XX63. ConservPredTFBS
		XX64. HumanEnhancer
		XX65. RNAedit
		XX66. PolyPhen2
		XX67. SIFT
		XX68. LS-SNP
		XX69. UniProt
		XX70. EqtlMethMetabStudy
		71. EFO string

	'''
	items = line.rstrip().split('\t')
	if len(items) < 61:
		assert False, line
	if items[11] in diseases or items[60] in efos:
		return GWAS_Association(
			pvalue = float(items[10]),
			snp = "rs" + items[4],
			disease = Disease(name = items[11], efo = items[70]),
			source = "GRASP",
			study = items[7]
		)
	else:
		return None

def Phewas_Catalog(diseases, efos):
	"""

		Returns all GWAS SNPs associated to a disease in PhewasCatalog
		Args:
		* [ string ] (trait descriptions)
		* [ string ] (trait EFO identifiers)
		Returntype: [ GWAS_Association ]

	"""
	file = open(DATABASES_DIR+"/Phewas_Catalog.txt")
	res = [ get_phewas_catalog_association(line, diseases, efos) for line in file ]
	res = filter(lambda X: X is not None, res)

	if DEBUG:
		print "\tFound %i GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in Phewas Catalog" % (len(res), ", ".join(diseases), ", ".join(efos))

	return res

def get_phewas_catalog_association(line, diseases, efos):
	'''

		Phewas Catalog format:
		1. chromosome
		2. snp
		3. phewas phenotype
		4. cases
		5. p-value
		6. odds-ratio
		7. gene_name
		8. phewas code
		9. gwas-associations
		10. [Inserte] EFO identifier (or N/A)

	'''
	items = line.rstrip().split('\t')
	if items[2] in diseases or items[9] in efos:
		return GWAS_Association (
			pvalue = float(items[4]),
			snp = items[1],
			disease = Disease(name = items[3], efo = items[9]), 
			source = "Phewas Catalog",
			study = None
		)
	else:
		return None

def GWAS_DB(diseases, efos):
	"""

		Returns all GWAS SNPs associated to a disease in GWAS_DB
		Args:
		* [ string ] (trait descriptions)
		* [ string ] (trait EFO identifiers)
		Returntype: [ GWAS_Association ]

	"""
	file = open(DATABASES_DIR+"/GWAS_DB.txt")

	# Note the format of the EFO strings is modified in this file format, so we need to change the queries
	efos2 = [re.sub("_", "ID:", efo) for efo in efos]
	res = [ get_gwas_db_association(line, diseases, efos) for line in file ]
	res = filter(lambda X: X is not None, res)

	if DEBUG:
		print "\tFound %i GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in GWAS DB" % (len(res), ", ".join(diseases), ", ".join(efos))

	return res

def get_gwas_db_association(line, diseases, efos):
	'''

		GWAS DB data
		1. CHR
		2. POS
		3. SNPID
		4. P_VALUE
		5. PUBMED ID
		6. MESH_TERM
		7. EFO_ID

	'''
	items = line.rstrip().split('\t')
	if items[5] in diseases or items[6] in efos:
		return GWAS_Association(
			pvalue = float(items[3]),
			snp = items[2],
			disease = Disease(name = items[5], efo = items[6]),
			source = "GWAS DB",
			study = items[4]
		)
	else:
		return None

def gwas_snps_to_genes(gwas_snps, populations, tissue_weights):
	"""

		Associates Genes to gwas_snps of interest
		Args:
		* [ GWAS_Association ]
		* { population_name: scalar (weight) }
		* { tissue_name: scalar (weight) }
		Returntype: [ GeneCluster_Association ]

	"""
	# Must set the tissue settings before separating out the gwas_snps
	if tissue_weights is None:
		tissue_weights = gwas_snps_to_tissue_weights(gwas_snps)

	clusters = cluster_gwas_snps(gwas_snps, populations)
	res = concatenate(cluster_to_genes(cluster, tissue_weights, populations) for cluster in clusters)

	if DEBUG:
		print "\tFound %i genes associated to all clusters" % (len(res))

	if len(res) == 0:
		return []
	else:
		# DEBUG
		# return [ sorted(res, key=lambda X: X.score)[-1] ]
		return sorted(res, key=lambda X: X.score)

def gwas_snps_to_tissue_weights(gwas_snps):
	"""

		Associates list of tissues to list of gwas_snps
		Args:
		* [ GWAS_SNP ]
		Returntype: [ string ]

	"""
	return ['Whole_Blood'] # See FORGE??

def cluster_gwas_snps(gwas_snps, populations):
	"""

		Bundle together gwas_snps within LD threshold
		* [ GWAS_SNP ]
		* [ Bio::Ensembl::Variation::Population ]
		Returntype: [ GWAS_Cluster ]

	"""
	gwas_snp_locations = get_gwas_snp_locations(gwas_snps)

	if DEBUG:
		"Found %i locations from %i GWAS SNPs" % (len(gwas_snp_locations), len(gwas_snps))

	preclusters = filter (lambda X: X is not None, [ gwas_snp_to_precluster(gwas_snp_location, populations) for gwas_snp_location in gwas_snp_locations ])
	clusters = merge_preclusters(preclusters)

	if DEBUG:
		"Found %i clusters from %i GWAS SNP locations" % (len(clusters), len(gwas_snp_locations))

	return clusters


def calculate_LD_window(snp, window_len=500000,population='CEPH',cutoff=0.5,db=0):

    """
    Given a SNP id, calculate the pairwise LD between all SNPs within window_size base pairs.
    """

    SNP_id = snp.rsID
    ### Get the SNP location from ENSEMBL
    position = int(snp.pos)

    ### Define the necessary region.
    from_pos = position - (window_len / 2)
    to_pos = position + (window_len / 2)
    chromosome = snp.chrom
    region = '{}:{}-{}'.format(chromosome,from_pos,to_pos)

    chrom_file = "ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf" % (chromosome)

    ### Extract this region out from the 1000 genomes BCF

    extract_region_comm = "bcftools view -r %s %s -O v -o region.vcf" % (region, os.path.join(DATABASES_DIR, '1000Genomes', population, chrom_file))

    subprocess.call(extract_region_comm.split(" "))
    region_file = open('region.vcf','r')
    region_vcf = region_file.read()


    ### Find the order of SNPs in the VCF
    SNPs_order = re.findall('rs[0-9]+', region_vcf)


    ### Calculate the pairwise LD using plink2
    plinkcomm = "plink --vcf region.vcf --r2 --ld-snp {} --inter-chr --out LDwindow".format(SNP_id)
    print plinkcomm
    plinkcomm_list = plinkcomm.split(" ")

    try:
        subprocess.call(plinkcomm_list)
    except:
        raise Exception("Plink2 needs to be installed")
        sys.exit()

    ### Remove intermediate region VCF file






    ### Remove intermediate LD file


    f = open('LDwindow.ld','r')
    g = f.read().splitlines()

    r2_dict = {}
    for s in [l.split() for l in g][1:]:
        ld = s[-1]
        the_snp = SNP(s[-2].split(';')[0], s[-4], int(s[-3]))
        r2_dict[the_snp] = float(ld)


    pruned_ld_snps = dict([x for x in r2_dict.items() if x[1] > cutoff])

    if db != 1:
        subprocess.call(['rm', 'region.vcf'])
    if db != 1:
        subprocess.call(['rm', 'LDwindow.ld', 'LDwindow.log', 'LDwindow.nosex', 'out.log'])

    return pruned_ld_snps


def gwas_snp_to_precluster(gwas_snp, populations):
    """

        Extract neighbourhood of GWAS snp
        Args:
        * [ GWAS_SNP ]
        * [ string ] (populations)
        Returntype: GWAS_Cluster

    """
    # Get all LD values around SNP
    snp = gwas_snp.snp
    mapped_ld_snps_dict = calculate_LD_window(snp)
    mapped_ld_snps = mapped_ld_snps_dict.keys()
    return GWAS_Cluster(gwas_snps = [ gwas_snp ],ld_snps = mapped_ld_snps)

def get_gwas_snp_locations(gwas_snps):
	"""

		Extract locations of GWAS SNP:
		Args
		* GWAS_SNP
		Returntype: [ GWAS_SNP ]

	"""
	original_gwas_snp = dict((gwas_snp.snp, gwas_snp) for gwas_snp in gwas_snps)
	mapped_snps = get_snp_locations(original_gwas_snp.keys())
	return [
		GWAS_SNP(
			snp = mapped_snp,
			pvalue = original_gwas_snp[mapped_snp.rsID].pvalue,
			evidence = original_gwas_snp[mapped_snp.rsID].evidence
		)
		for mapped_snp in mapped_snps
	]

def merge_preclusters(preclusters):
	"""

		Bundle together preclusters that share one LD snp
		* [ Cluster ]
		Returntype: [ Cluster ]

	"""
	snp_owner = dict()
	kill_list = list()
	for cluster in preclusters:
		for ld_snp in cluster.ld_snps:
			if ld_snp in snp_owner and snp_owner[ld_snp] is not cluster:
				other_cluster = snp_owner[ld_snp]
				print "Overlap between %i and %i" % (id(cluster), id(other_cluster))

				# Merge data from current cluster into previous cluster
				merged_gwas_snps = other_cluster.gwas_snps + cluster.gwas_snps
				merged_ld_snps = dict((ld_snp.rsID, ld_snp) for ld_snp in cluster.ld_snps + other_cluster.ld_snps).values()
				snp_owner[ld_snp] = GWAS_Cluster(
					gwas_snps = merged_gwas_snps,
					ld_snps = merged_ld_snps
				)

				# Mark for deletion
				kill_list.append(cluster)
				for snp in cluster.ld_snps:
					snp_owner[snp] = other_cluster

				# Exit from that cluster
				break
			else:
				snp_owner[ld_snp] = cluster

	res = filter(lambda cluster: cluster not in kill_list, preclusters)

	if DEBUG:
		print "\tFound %i clusters from the GWAS peaks" % (len(res))

	return res


def cluster_to_genes(cluster, tissues, populations):
    """

        Associated Genes to a cluster of gwas_snps
        Args:
        * [ Cluster ]
        * { tissue_name: scalar (weights) }
        * { population_name: scalar (weight) }
        Returntype: [ GeneCluster_Association ]

    """
    # Obtain interaction data from LD snps
    associations = ld_snps_to_genes(cluster.ld_snps, tissues)

    # Compute LD based scores
    top_gwas_hit = sorted(cluster.gwas_snps, key=lambda X: X.pvalue)[-1]
    ld = get_lds_from_top_gwas(top_gwas_hit.snp, cluster.ld_snps, populations)
    pics = PICS(ld, top_gwas_hit.pvalue)

    gene_scores = dict(
        ((association.gene, association.snp), (association, association.score * ld[association.snp]))
        for association in associations if association.snp in ld
    )

    if len(gene_scores) == 0:
        return []

    # OMIM exception
    max_score = max(X[1] for X in gene_scores.values())
    for gene, snp in gene_scores:
        if len(gene_to_phenotypes(gene)):
            gene_scores[(gene, snp)][1] = max_score


    res = [
        GeneCluster_Association(
            gene = gene,
            score = total_score(pics[snp], gene_scores[(gene, snp)][1]),
            cluster = cluster,
            evidence = gene_scores[(gene, snp)][:1] # This is a [ GeneSNP_Association ]
        )
        for (gene, snp) in gene_scores if snp in pics
    ]

    if DEBUG:
        print "\tFound %i genes associated around GWAS SNP %s" % (len(res), top_gwas_hit.snp.rsID)

    # Pick the association with the highest score
    # DEBUG
    #return [ sorted(res, key=lambda X: X.score)[-1] ]
    return sorted(res, key=lambda X: X.score)

def get_lds_from_top_gwas(gwas_snp, ld_snps, population='CEPH', region=None,db=0, cutoff=0.5):
    """
    For large numbers of SNPs, best to specify SNP region with chrom:to-from, e.g. 1:7654947-8155562
    For small numbers (<10), regions are extracted from ENSEMBL REST API.
    SNPs can be inputted in a list or from a file with one SNP id per line.
    """


    if gwas_snp.rsID not in [x.rsID for x in ld_snps]:
        ld_snps.append(gwas_snp)


    snp_position_map = {}
    for s in ld_snps:
        snp_position_map[s.rsID] = s.pos



    assert [x.chrom for x in ld_snps] == ([ld_snps[0].chrom] * len(ld_snps))
    positions = [x.pos for x in ld_snps]
    region = '{}:{}-{}'.format(ld_snps[0].chrom, min(positions), max(positions))

    chromosome = ld_snps[0].chrom


    SNPs_filepath = 'snp_rsIDs.txt'
    h = open('snp_rsIDs.txt', 'w')
    for x in ld_snps:
        h.write('{}\n'.format(str(x.rsID)))
    h.close()


    ### Extract the required region from the VCF
    chrom_file = 'ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf' % (chromosome)

    extract_region_comm = "bcftools view -r %s %s -O z -o region.vcf.gz" % (region, os.path.join(DATABASES_DIR, '1000Genomes', population, chrom_file))
    subprocess.call(extract_region_comm.split(" "))
    region_file = "region.vcf.gz"


    ### Extract the list of SNP ids from this region
    vcfcomm = "vcftools --gzvcf {} --snps {} --recode --stdout".format(region_file, SNPs_filepath)
    print vcfcomm
    vcf = subprocess.check_output(vcfcomm.split(" "))


    ### Remove intermediate region VCF file


    f = open('snps.vcf', 'w')
    f.write(vcf)
    f.close()

    ### Extract out the order of SNPs
    SNPs_order = re.findall('rs[0-9]+', vcf)


    ### Use plink2 to calculate pairwise LD between these SNPs.
    plinkcomm = "plink --vcf snps.vcf --r2 square --out LD"
    print plinkcomm
    plinkcomm_list = plinkcomm.split(" ")
    subprocess.call(plinkcomm_list)



    ### Read from the generated results file and output an array.
    LD_file = open('LD.ld','r')
    g = LD_file.read()
    LD_array = [x.split('\t') for x in g.splitlines()]
    LD_file.close


    snp_index = SNPs_order.index(gwas_snp.rsID)
    ld_vector = LD_array[snp_index]


    SNPs = [SNP(snp_id, chromosome, snp_position_map[snp_id]) for snp_id in SNPs_order]
    ld_scores = [(SNPs[i],float(ld_vector[i])) for i in range(len(ld_vector))]
    ld_scores_cutoff = dict([x for x in ld_scores if x[1] > cutoff])

    ### Remove intermediate files
    if db != 1:
        subprocess.call(['rm', 'LD.ld', 'LD.log', 'LD.nosex', 'out.log'])
        subprocess.call(['rm', 'snps.vcf'])
        subprocess.call(['rm', 'region.vcf.gz'])
        subprocess.call(['rm', 'snp_rsIDs.txt'])

    return ld_scores_cutoff



def PICS(ld, pvalue):
    minus_log_pvalue = - math.log(pvalue) / math.log(10);
    SD = dict()
    Mean = dict()
    prob = dict()
    sum = 0

    for snp in ld.keys():
        if snp in ld:
            SD[snp] = math.sqrt(1 - math.sqrt(ld[snp]) ** 6.4) * math.sqrt(minus_log_pvalue) / 2;
            Mean[snp] = ld[snp] * minus_log_pvalue;
        else:
            # Defaults for remote SNPs
            SD[snp] = 0
            Mean[snp] = 1 + minus_log_pvalue

        # Calculate the probability of each SNP
        if SD[snp]:
            prob[snp] = 1 - pnorm(minus_log_pvalue, Mean[snp], SD[snp])
        else:
            prob[snp] = 1

        # Normalisation sum
        sum += prob[snp]

    # Normalize the probabilies so that their sum is 1.
    return dict((snp, prob[snp] / sum) for snp in prob.keys())

def pnorm(x, mu, sd):
	"""

		Normal distribution PDF
		Args:
		* scalar: variable
		* scalar: mean
		* scalar: standard deviation
		Return type: scalar (probability density)

	"""
	return math.exp ( - ((x - mu) / sd) ** 2 / 2 ) / (sd * 2.5)

def total_score(pics, gene_score):
	"""

		Computes a weird mean function from ld_snp PICs score and Gene/SNP association score
		Args:
		* PICS: scalar
		* gene_score: scalar
		Returntype: scalar

	"""
	A = pics * (pics ** (1/3))
	B = gene_score * (gene_score ** (1/3))
	return ((A + B) / 2) ** 3

def rsIDs_to_genes(snp, tissues):
	"""

		Associates genes to single SNP
		Args:
		* SNP
		* [ string ] (tissues)
		Returntype: [ GeneSNP_Association ]

	"""
	if tissues is None:
		tissues = ["Whole_Blood"]
	return ld_snps_to_genes(get_snp_locations([snp]), tissues)

def ld_snps_to_genes(ld_snps, tissues):
	"""

		Associates genes to LD linked SNP
		Args:
		* [ SNP ]
		* [ string ] (tissues)
		Returntype: [ GeneSNP_Association ]

	"""
	# Search for SNP-Gene pairs:
	cisreg = cisregulatory_evidence(ld_snps, tissues)

	# Extract list of relevant SNPs:
	selected_snps = set(gene_snp_pair[1] for gene_snp_pair in cisreg)

	# Fallback solution: add nearest gene for dissappointing SNPs:
	if len(selected_snps) != len(ld_snps):
		nearest_gene_hash = nearest_genes(ld_snps)
		for snp in ld_snps:
			if snp not in selected_snps and snp in nearest_gene_hash:
				gene = nearest_gene_hash[snp.rsID]
				cisreg[(gene, snp)].append(Cisregulatory_Evidence(
						gene = gene,
						snp = snp,
						score = 1,
						source = 'Nearest',
						study = None,
						tissue = None))

	if len(selected_snps) > 0:
		# Extract SNP specific info:
		reg = regulatory_evidence(selected_snps, tissues)
	else:
		reg = []

	return [
		GeneSNP_Association(
			gene = gene,
			snp = snp,
			cisregulatory_evidence = cisreg[(gene, snp)],
			regulatory_evidence = reg[snp],
			score = sum(float(evidence.score) for evidence in reg[snp] + cisreg[(gene, snp)])
		)
		for (gene, snp) in cisreg
	]

def cisregulatory_evidence(ld_snps, tissues):
	"""

		Associates genes to LD linked SNP
		Args:
		* [ SNP ]
		* [ string ] (tissues)
		Returntype: { (Gene, SNP): Cisregulatory_Evidence }

	"""
	if DEBUG:
		print "Searching for cis-regulatory data on %i SNPs in all databases" % (len(ld_snps))
	evidence = concatenate(function(ld_snps, tissues) for function in ld_snp_to_gene_functions)

	filtered_evidence = filter(lambda association: association.gene.biotype == "protein_coding", evidence)

	# Group by (gene,snp) pair:
	res = collections.defaultdict(list)
	for association in filtered_evidence:
		res[(association.gene, association.snp)].append(association)

	if DEBUG:
		print "Found %i cis-regulatory interactions in all databases" % (len(res))

	return res

def GTEx(ld_snps, tissues):
	"""

		Returns all genes associated to a set of SNPs in GTEx
		Args:
		* [ SNP ]
		* [ string ] (tissues)
		Returntype: [ Cisregulatory_Evidence ]

	"""
	# Find all genes with 1Mb
	start = min(snp.pos for snp in ld_snps)
	end = max(snp.pos for snp in ld_snps)
	chrom = ld_snps[0].chrom

	server = 'http://grch37.rest.ensembl.org'
	ext = '/overlap/region/%s/%s:%i-%i?feature=gene;content-type=application/json' % (SPECIES, chrom, max(0, start - 1e6), end + 1e6)
	genes = [ Gene(
			name = gene['external_name'],
			id = gene['id'],
			chrom = gene['seq_region_name'],
			tss = int(gene['start']) if gene['strand'] > 0 else int(gene['end']),
			biotype = gene['biotype']
		)
		for gene in get_rest_json(server, ext)
	]
	snp_hash = dict( (snp.rsID, snp) for snp in ld_snps)
	res = concatenate((GTEx_gene(gene, tissues, snp_hash) for gene in genes))

	if DEBUG:
		print "\tFound %i interactions in GTEx" % (len(res))

	return res

def GTEx_gene(gene, tissues, snp_hash):
	"""

		Returns all SNPs associated to a gene in GTEx
		Args:
		* [ SNP ]
		* [ string ] (tissues)
		* { rsID: rsID }
		Returntype: [ Cisregulatory_Evidence ]

	"""
	res = concatenate(GTEx_gene_tissue(gene, tissue, snp_hash) for tissue in tissues)

	if DEBUG:
		print "\tFound %i genes associated to gene %s in GTEx" % (len(res), gene.id)

	return res

def GTEx_gene_tissue(gene, tissue, snp_hash):
	"""

		Returns all SNPs associated to a gene in GTEx in a given tissue
		Args:
		* [ SNP ]
		* [ string ] (tissues)
		* { rsID: rsID }
		Returntype: [ Cisregulatory_Evidence ]

	"""


	server = "http://rest.ensembl.org"
	ext = "/eqtl/id/%s/%s?content-type=application/json;statistic=p-value;tissue=%s" % ('homo_sapiens', gene.id, tissue);
	try:
		eQTLs = get_rest_json(server, ext)

		'''
			Example return object:
			[
				{
					'value': '0.804108648395327',
					'snp': 'rs142557973'
				},
			]
		'''

		res = [
			Cisregulatory_Evidence(
				snp = snp_hash[eQTL['snp']],
				gene = gene,
				tissue = tissue,
				score = 1,
				source = "GTEx",
				study = None
			)
			for eQTL in eQTLs if eQTL['snp'] in snp_hash
			if eQTL['value'] < 2.5e-5 
		]

		if DEBUG:
			print "\tFound %i SNPs associated to gene %s in tissue %s in GTEx" % (len(res), gene.id, tissue)

		return res
	except:
		return None

def chunks(l, n):
	for i in range(0, len(l), n):
		yield l[i:i+n]

def VEP(ld_snps, tissues):
	"""

		Returns all genes associated to a set of SNPs in VEP
		Args:
		* [ SNP ]
		* [ string ] (tissues)
		Returntype: [ Regulatory_Evidence ]

	"""
	server = "http://grch37.rest.ensembl.org"
	ext = "/vep/%s/id" % (SPECIES)
	list = concatenate(get_rest_json(server, ext, data = {"ids" : [snp.rsID for snp in chunk]}) for chunk in chunks(ld_snps, 999))
	'''

		Example output from VEP:
		[
			{
				'colocated_variants': [
					{
						'phenotype_or_disease': 1,
						'seq_region_name': '9',
						'eas_allele': 'C',
						'amr_maf': '0.4553',
						'strand': 1,
						'sas_allele': 'C',
						'id': 'rs1333049',
						'allele_string': 'G/C',
						'sas_maf': '0.4908',
						'amr_allele': 'C',
						'minor_allele_freq': '0.4181',
						'afr_allele': 'C',
						'eas_maf': '0.5367',
						'afr_maf': '0.2133',
						'end': 22125504,
						'eur_maf': '0.4722',
						'eur_allele': 'C',
						'minor_allele': 'C',
						'pubmed': [
							24262325,
						],
						'start': 22125504
					}
				],
				'assembly_name': 'GRCh38',
				'end': 22125504,
				'seq_region_name': '9',
				'transcript_consequences': [
					{
						'gene_id': 'ENSG00000240498',
						'variant_allele': 'C',
						'distance': 4932,
						'biotype': 'antisense',
						'gene_symbol_source': 'HGNC',
						'consequence_terms': [
							'downstream_gene_variant'
						],
						'strand': 1,
						'hgnc_id': 'HGNC:34341',
						'gene_symbol': 'CDKN2B-AS1',
						'transcript_id': 'ENST00000584020',
						'impact': 'MODIFIER'
					},
				],
				'strand': 1,
				'id': 'rs1333049',
				'most_severe_consequence': 'downstream_gene_variant',
				'allele_string': 'G/C',
				'start': 22125504
			}
		]

	'''

	snp_hash = dict( (snp.rsID, snp) for snp in ld_snps)
	transcript_consequences = filter(lambda X: 'transcript_consequences' in X, list)
	res = [
		Cisregulatory_Evidence(
			snp = snp_hash[hit['input']],
			gene = get_ensembl_gene(consequence['gene_id']),
			score = VEP_impact_to_score[consequence['impact']],
			source = "VEP",
			study = None,
			tissue = None
		)
		for hit in transcript_consequences for consequence in hit['transcript_consequences']
	]

	if DEBUG:
		print "\tFound %i interactions in VEP" % (len(res))

	return res

def Fantom5(ld_snps, tissues):
	"""

		Returns all genes associated to a set of SNPs in Fantom5
		Args:
		* [ SNP ]
		* [ string ] (tissues)
		Returntype: [ Regulatory_Evidence ]

	"""
	intersection = overlap_snps_to_bed(ld_snps, DATABASES_DIR + "/Fantom5.bed")
	fdr_model = pickle.load(open(DATABASES_DIR + "/Fantom5.fdrs"))
	snp_hash = dict( (snp.rsID, snp) for snp in ld_snps)
	hits  = filter(lambda X: X is not None, map(lambda X: get_fantom5_evidence(X, fdr_model, snp_hash), intersection))

	if DEBUG:
		print "\tFound %i overlaps with %i SNPs to Fantom5" % (len(hits), len(ld_snps))

	res = filter(lambda X: X.score, hits)

	if DEBUG:
		print "\tFound %i interactions in Fantom5" % (len(res))

	return res

def get_fantom5_evidence(feature, fdr_model, snp_hash):
	'''
		Parse output: first 12 columns are from Fantom5 file, the next 4 are LD SNP coords
		1. chrom
		2. chromStart
		3. chromEnd
		4. HGNC
		5. FDR:fdr

		6. chrom
		7. start
		8. end
		9. rsID
	'''
	gene = get_gene(feature[3])
	if gene is None:
		return None
	snp = snp_hash[feature[8]]
	score = STOPGAP_FDR(snp, gene, fdr_model)

	return Cisregulatory_Evidence(
			snp = snp,
			gene = gene,
			source = "Fantom5",
			score = score,
			study = None,
			tissue = None
		)

def DHS(ld_snps, tissues):
	"""

		Returns all genes associated to a set of SNPs in DHS
		Args:
		* [ SNP ]
		* [ string ] (tissues)
		Returntype: [ Regulatory_Evidence ]

	"""
	intersection = overlap_snps_to_bed(ld_snps, DATABASES_DIR + "/DHS.bed")
	fdr_model = pickle.load(open(DATABASES_DIR+"/DHS.fdrs"))
	snp_hash = dict( (snp.rsID, snp) for snp in ld_snps)
	res = filter (lambda X: X is not None and X.score, (get_dhs_evidence(feature, fdr_model, snp_hash) for feature in intersection))

	if DEBUG:
		print "\tFound %i gene associations in DHS" % len(res)

	return res

def get_dhs_evidence(feature, fdr_model, snp_hash):
	"""

		Parse output of bedtools searching for DHS evidence
		Args:
		* string
		Returntype: Regulatory_Evidence

	"""
	'''

		First 5 columns are from DHS file, the next 4 are LD SNP coords
		Format of DHS correlation files:
		1. chrom
		2. chromStart
		3. chromEnd
		4. HGNC
		5. Correlation

		6. chrom
		7. start
		8. end
		9. rsID

	'''
	gene = get_gene(feature[3])
	if gene is None:
		return None
	snp = snp_hash[feature[8]]
	score = STOPGAP_FDR(snp, gene, fdr_model)

	return Cisregulatory_Evidence(
		snp = snp,
		gene = gene,
		score = score,
		source = "DHS",
		study = None,
		tissue = None
	)

def PCHIC(ld_snps, tissues):
	"""

		Returns all genes associated to a set of SNPs in PCHIC
		Args:
		* [ SNP ]
		* [ string ] (tissues)
		Returntype: [ Regulatory_Evidence ]

	"""
	intersection = overlap_snps_to_bed(ld_snps, DATABASES_DIR + "/pchic.bed")
	snp_hash = dict( (snp.rsID, snp) for snp in ld_snps)
	res = filter (lambda X: X is not None and X.score, (get_pchic_evidence(feature, snp_hash) for feature in intersection))

	if DEBUG:
		print "\tFound %i gene associations in PCHIC" % len(res)

	return res

def get_pchic_evidence(feature, snp_hash):
	"""

		Parse output of bedtools searching for PCHIC evidence
		Args:
		* string
		Returntype: Regulatory_Evidence

	"""
	'''

		First 5 columns are from PCHIC file, the next 4 are LD SNP coords
		Format of PCHIC correlation files:
		1. chrom
		2. chromStart
		3. chromEnd
		4. score
		5. gene_id

		6. chrom
		7. start
		8. end
		9. rsID

	'''
	gene = get_gene(feature[4])
	if gene is None:
		return None
	snp = snp_hash[feature[8]]

	return Cisregulatory_Evidence(
		snp = snp,
		gene = gene,
		score = 1,
		source = "PCHIC",
		study = None,
		tissue = None
	)

def nearest_genes(snps):
	"""

		Return nearest gene to SNP
		Args:
		* [ SNP ]
		Returntype: dict(rsID => Gene)

	"""
	snps = list(snps)
	bed = DATABASES_DIR + "/Ensembl_TSSs.bed"
	SNP_string = "\n".join("\t".join((snp.chrom, str(snp.pos), str(snp.pos+1), snp.rsID)) for snp in snps)
	SNP_bt = pybedtools.BedTool(SNP_string, from_string=True)

	if not os.path.isfile(bed + '.gz') or not os.path.isfile(bed + '.gz.tbi'):
		Annotation_bt = pybedtools.BedTool(bed)
		Annotation_bt_indexed = Annotation_bt.tabix()
	else:
		Annotation_bt_indexed = pybedtools.BedTool(bed + '.gz')

	'''
		Parse output: first 4 columns are from SNP file, next 4 Gene coords
		1. chrom
		2. start
		3. end
		4. rsID

		5. chrom
		6. chromStart
		7. chromEnd
		8. HGNC

	'''
	return dict((row[3], get_gene(row[7])) for row in SNP_bt.closest(Annotation_bt_indexed, wa=True, wb=True))

def STOPGAP_FDR(snp, gene, fdr_model):
	"""

		Special function for cis-regulatory interactions
		Args:
		* rsID
		* ENSG stable ID
		*
		Returntype: scalar

	"""

	if gene.chrom != snp.chrom:
		return 0

	distance = abs(snp.pos - gene.tss)

	if distance > fdr_model.MAX_DISTANCE:
		return 0

	bin = int(distance / fdr_model.BIN_WIDTH)

	if bin not in fdr_model.FDR:
		return 0

	FDR = fdr_model.FDR[bin]

	if FDR is None:
		return 0
	elif FDR < .6:
		return 2
	elif FDR < .85:
		return 1
	else:
		return 0

def get_gene(gene_name):
	"""

		Get gene details from name
		* string
		Returntype: Gene

	"""
	if gene_name not in known_genes:
		if gene_name[:4] != 'ENSG':
			known_genes[gene_name] = fetch_gene(gene_name)
		else:
			known_genes[gene_name] = fetch_gene_id(gene_name)
	return known_genes[gene_name]

def fetch_gene(gene_name):
	"""

		Get gene details from name
		* string
		Returntype: Gene

	"""
	server = "http://grch37.rest.ensembl.org"
	ext = "/lookup/symbol/%s/%s?content-type=application/json;expand=1" % (SPECIES, gene_name)
	try:
		hash = get_rest_json(server, ext)
		return Gene(
			name = gene_name,
			id = hash['id'],
			chrom = hash['seq_region_name'],
			tss = int(hash['start']) if hash['strand'] > 0 else int(hash['end']),
			biotype = hash['biotype']
			)
	except:
		return None

def fetch_gene_id(gene_id):
	"""

		Get gene details from name
		* string
		Returntype: Gene

	"""
	server = "http://grch37.rest.ensembl.org"
	ext = "/lookup/id/%s?content-type=application/json;expand=1" % (gene_id)
	try:
		hash = get_rest_json(server, ext)
		return Gene(
			name = hash['display_name'],
			id = hash['id'],
			chrom = hash['seq_region_name'],
			tss = int(hash['start']) if hash['strand'] > 0 else int(hash['end']),
			biotype = hash['biotype']
			)
	except:
		return None


def get_ensembl_gene(ensembl_id):
	"""

		Get gene details from name
		* string
		Returntype: Gene

	"""
	if ensembl_id not in known_genes:
		known_genes[ensembl_id] = fetch_ensembl_gene(ensembl_id)
	return known_genes[ensembl_id]

def fetch_ensembl_gene(ensembl_id):
	"""

		Get gene details from name
		* string
		Returntype: Gene

	"""
	server = "http://grch37.rest.ensembl.org"
	ext = "/lookup/id/%s?content-type=application/json;expand=1" % (ensembl_id)
	hash = get_rest_json(server, ext)
	return Gene(
		name = hash['display_name'],
		id = ensembl_id,
		chrom = hash['seq_region_name'],
		tss = int(hash['start']) if hash['strand'] > 0 else int(hash['end']),
		biotype = hash['biotype']
		)

def get_snp_locations(rsIDs):
	"""

		Get SNP details from rsID
		* string
		Returntype: SNP

	"""

	server = "http://grch37.rest.ensembl.org"
	ext = "/variation/%s?content-type=application/json" % (SPECIES)
	hash = concatenate_hashes(get_rest_json(server, ext, data={'ids':chunk}) for chunk in chunks(rsIDs, 999))

	'''
		Example response:
		{
			"rs56116432": {
				"source": "Variants (including SNPs and indels) imported from dbSNP",
				"mappings": [
					{
						"location": "9:133256042-133256042",
						"assembly_name": "GRCh38",
						"end": 133256042,
						"seq_region_name": "9",
						"strand": 1,
						"coord_system": "chromosome",
						"allele_string": "C/T",
						"start": 133256042
					},
				],
				"name": "rs56116432",
				"MAF": "0.00259585",
				"ambiguity": "Y",
				"var_class": "SNP",
				"synonyms": [],
				"evidence": [
					"Multiple_observations",
					"Frequency",
					"1000Genomes",
					"ESP",
					"ExAC"
				],
				"ancestral_allele": "C",
				"minor_allele": "T",
				"most_severe_consequence": "missense_variant"
			}
		}
	'''
	return [
		SNP(
			rsID = rsID,
			chrom = mapping['seq_region_name'],
			pos = (int(mapping['start']) + int(mapping['end'])) / 2
		)
		for rsID in hash
		for mapping in hash[rsID]['mappings']
	]

def regulatory_evidence(snps, tissues):
	"""

		Extract regulatory evidence linked to SNPs and stores them in a hash
		* [ SNP ]
		* [ string ]
		Returntype: [ Regulatory_Evidence ]

	"""
	if DEBUG:
		print "Searching for regulatory data on %i SNPs in all databases" % (len(snps))

	res = concatenate(function(snps, tissues) for function in snp_regulatory_functions)

	if DEBUG:
		print "Found %i regulatory SNPs among %i in all databases" % (len(res), len(snps))

	# Group by SNP
	hash = collections.defaultdict(list)
	for hit in res:
		hash[hit.snp].append(hit)

	return hash

def GERP(snps, tissues):
	"""

		Extract GERP score at position
		Args:
		* [ SNP ]
		Returntype: [ Regulatory_Evidence ]

	"""
	return map(GERP_at_snp, snps)

def GERP_at_snp(snp):
	"""

		Extract GERP score at position
		Args:
		* rsID
		Returntype: Regulatory_Evidence
	"""
	server = "http://grch37.rest.ensembl.org"
	ext = "/vep/human/id/%s?content-type=application/json;Conservation=1" % snp.rsID
	obj = get_rest_json(server, ext)
	return Regulatory_Evidence(
			snp = snp,
			score = float(obj['Conservation']), # TODO: score on FDR?
			source = 'GERP'
	)

def Regulome(ld_snps, tissues):
	"""

		Extract Regulome score at sns of interest
		Args:
		* [ SNP ]
		* [ string ]
		Returntype: [ Regulatory_Evidence ]

	"""
	snp_hash = dict( (snp.rsID, snp) for snp in ld_snps)
	intersection = overlap_snps_to_bed(ld_snps, DATABASES_DIR + "/Regulome.bed")
	res = filter (lambda X: X.score, (get_regulome_evidence(feature, snp_hash) for feature in intersection))

	if DEBUG:
		print "\tFound %i regulatory variants in Regulome" % (len(res))

	return res

def overlap_snps_to_bed(ld_snps, bed):
    '''
        Find overlaps between SNP elements and annotated Bed file
        Args:
        * [ SNP ]
        * string (location of bed file)
        Returntype: pybedtools Interval iterator

    '''
    ld_snps = list(ld_snps)
    SNP_string = "\n".join(( "\t".join((snp.chrom, str(snp.pos), str(snp.pos+1), snp.rsID)) for snp in ld_snps ))
    SNP_bt = pybedtools.BedTool(SNP_string, from_string=True)
    max_pos = max([x.pos for x in ld_snps])
    min_pos = min([x.pos for x in ld_snps])
    chrom = ld_snps[0].chrom

    if not os.path.isfile(bed + '.gz') or not os.path.isfile(bed + '.gz.tbi'):
        Annotation_bt = pybedtools.BedTool(bed)
        Annotation_bt_indexed = Annotation_bt.tabix()
    else:
        Annotation_bt_indexed = pybedtools.BedTool(bed + '.gz')

    intersection = Annotation_bt_indexed.tabix_intervals('{}:{}-{}'.format(chrom,min_pos,max_pos)).intersect(SNP_bt, wa=True, wb=True)
    return intersection

def get_regulome_evidence(feature, snp_hash):
	"""

		Extract Regulome score from bedtools output
		Args:
		* string
		Returntype: Regulatory_Evidence

	"""

	'''
	First 4 columns are from Regulome file, the next 4 are LD SNP coords
	Regulome file format:
	1. chrom
	2. start
	3. end
	4. category

	LD SNP coords:
	5. chrom
	6. start
	7. end
	8. rsID
	'''
	return Regulatory_Evidence(
			snp = snp_hash[feature[7]],
			source = "Regulome",
			score = 2 if feature[3][0] == '1' or feature[3][0] == '2' else 1,
			study = None,
			tissue = None
		)

def gene_to_MeSH(gene):
	"""

		Look up MeSH annotations for gene
		Args:
		* [ string ] (gene names)
		Return type: [ string ] (annotations)
	"""
	server = "http://gene2mesh.ncibi.org"
	ext = "/fetch?genesymbol=%s&taxid=%s" % (gene.name, NCBI_Taxon_ID[SPECIES])
	print '>>>>>>>>>>>>>>>>>>'
	print str(server)+str(ext)
	print '>>>>>>>>>>>>>>>>>>'
	response = requests.get(str(server)+str(ext))

	if not response.ok:
		sys.stderr.write("Failed to get proper response to query %s%s\n" % (server, ext) )
		sys.stderr.write(response.content + "\n")
		response.raise_for_status()
		print repr(response)

	'''
	Example MeSH output:

	{
		'Gene2MeSH': {
			'Request': {
				'ParameterSet': {
					'Tool': 'none',
					'GeneSymbol': 'csf1r',
					'TaxonomyID': '9606',
					'Limit': '1000',
					'Email': 'anonymous'
				},
				'type': 'fetch'
			},
			'Response': {
				'ResultSet': {
					'Result': [
						{
							'FisherExact': {
								'content': '1.8531319238671E-230',
								'type': 'p-value'
							},
							'ChiSquare': '112213.6506462',
							'Fover': '1498.1813411401',
							'MeSH': {
								'Qualifier': {
									'Name': 'metabolism'
								},
								'Descriptor': {
									'TreeNumber': [
										'D08.811.913.696.620.682.725.400.500',
										'D12.776.543.750.060.492',
										'D12.776.543.750.705.852.150.150',
										'D12.776.543.750.750.400.200.200',
										'D12.776.624.664.700.800',
										'D23.050.301.264.035.597',
										'D23.101.100.110.597'
									],
									'Identifier': 'D016186',
									'Name': 'Receptor, Macrophage Colony-Stimulating Factor',
									'UMLSID': {}
								}
							},
							'DocumentSet': {
								'type': 'pubmed',
								'PMID': [
								]
							},
							'Gene': {
								'Taxonomy': {
									'Identifier': '9606'
								},
								'Identifier': '1436',
								'type': 'Entrez',
								'Description': 'colony stimulating factor 1 receptor',
								'Symbol': 'CSF1R'
							}
						},
					],
					'sort': 'FisherExact',
					'count': '94',
					'order': 'ascending'
				},
				'Copyright': {
					'Details': 'http://nlp.ncibi.org/Copyright.txt',
					'Year': '2009',
					'Statement': 'Copyright 2009 by the Regents of the University of Michigan'
				},
				'Support': {
					'Details': 'http://www.ncibi.org',
					'GrantNumber': 'U54 DA021519',
					'Statement': 'Supported by the National Institutes of Health as part of the NIH\\\'s National Center for Integrative Biomedical Informatics (NCIBI)'
				}
			}
		}
	}
	'''

	print response.content
	if len(response.content):
		try:
			hash = xmltodict.parse(response.content)
			print repr(hash)
			hits = hash['Gene2MeSH']['Request']['ResultSet']['Result']
			# XML subtletly: if a tag is repeated, a list of objects is produced,
			# else a single object. Careful when unpacking!
			if hits is list:
				return [hit['MeSH']['Descriptor']['Name'] for hit in hits]
			else:
				return [hits['MeSH']['Descriptor']['Name']]
		except:
			return []
	else:
		return []

def gene_to_phenotypes(gene):
	"""

		Look up phenotype annotations for gene
		Args:
		* ENSG stable ID
		Return type: [ OMIM Phenotype ]

	"""
	return [] # TODO and remove stopper
	if gene not in phenotype_cache:
		phenotype_cache[gene['stable_id']] = get_gene_phenotypes(gene)

	return phenotype_cache[gene['stable_id']]

def get_rest_json(server, ext, data=None):
	"""
		Args:
		* String (server name)
		* String (extension string)
		Return type: JSON object

	"""
	retries = 0

	while True:
		if data is None:
			headers = { "Content-Type" : "application/json" }
			r = requests.get(str(server)+str(ext), headers = headers)
		else:
			headers = {'Content-Type': 'application/json', 'Accept': 'application/json'}
			r = requests.post(str(server)+str(ext), headers = headers, data = json.dumps(data))

		if DEBUG:
			sys.stderr.write("REST JSON Query: %s%s\n" % (server, ext))

		if not r.ok:
			sys.stderr.write("Failed to get proper response to query %s%s\n" % (server, ext) )
			sys.stderr.write("With headers:\n" + repr(headers) + "\n")
			if data is not None:
				sys.stderr.write("With data:\n" + repr(data) + "\n")
			if 'Retry-After' in r.headers:
				time.sleep(int(r.headers['Retry-After']))
				retries += 1
				continue
			r.raise_for_status()
			sys.exit()

		try:
			return r.json()
		except:
			sys.stderr.write("Failed to get proper response to query %s%s\n" % (server, ext) )
			raise

	# Failed too many times
	sys.exit()

def concatenate(list):
	"""

		Shorthand to concatenate a list of lists
		Args: [[]]
		Returntype: []

	"""
	return sum(filter(lambda elem: elem is not None, list), [])

def concatenate_hashes(list):
	"""

		Shorthand to concatenate a list of lists
		Args: [[]]
		Returntype: []

	"""
	return dict(sum(map(lambda X: X.items(), filter(lambda elem: elem is not None, list)), []))

def pretty_snp_output(associations):
	"""

		Prints association stats in roughly the same format as STOPGAP for a cluster of SNPs
		Args: [ GeneSNP_Association ]
		Returntype: String

	"""
	header = "\t".join(['snp_rsID', 'gene_symbol', 'gene_id', 'score', 'VEP', 'Regulome', 'GTEx', 'DHS', 'Fantom5', 'PCHIC', 'Nearest'])
	content = map(pretty_snp_association, associations)
	return "\n".join([header] + content) + "\n"

def pretty_snp_association(association):
	"""

		Prints association stats in roughly the same format as STOPGAP for a cluster of SNPs
		Args: GeneCluster_Association
		Returntype: String

	"""
	snp = association.snp
	gene_name = association.gene.name
	gene_id = association.gene.id
	score = association.score

	functional_scores = collections.defaultdict(int)
	for evidence in association.regulatory_evidence + association.cisregulatory_evidence:
		functional_scores[evidence.source] += evidence.score
	
	results = [snp.rsID, gene_name, gene_id, str(score)]
	results += [str(functional_scores[functional_source]) for functional_source in ['VEP', 'Regulome', 'GTEx', 'DHS', 'Fantom5', 'PCHIC', 'Nearest']]
	return "\t".join(results)


def pretty_output(associations):
	"""

		Prints association stats in roughly the same format as STOPGAP for a cluster of SNPs
		Args: [ GeneCluster_Association ]
		Returntype: String

	"""
	header = "\t".join(['ld_snp_rsID', 'gene_symbol', 'gene_id', 'disease_names', 'disease_efo_ids', 'score', 'gwas_snp_ids', 'GWAS Catalog', 'GRASP', 'GWAS DB', 'Phewas Catalog', 'VEP', 'Regulome', 'GTEx', 'DHS', 'Fantom5', 'PCHIC', 'Nearest'])
	content = map(pretty_cluster_association, associations)
	return "\n".join([header] + content)

def pretty_cluster_association(association):
	"""

		Prints association stats in roughly the same format as STOPGAP for a cluster of SNPs
		Args: GeneCluster_Association
		Returntype: String

	"""
	gene_name = association.gene.name
	gene_id = association.gene.id
	cluster = association.cluster
	gwas_snps = cluster.gwas_snps
	disease_names = list(set(gwas_association.disease.name for gwas_snp in gwas_snps for gwas_association in gwas_snp.evidence))
	disease_efos = list(set(gwas_association.disease.efo for gwas_snp in gwas_snps for gwas_association in gwas_snp.evidence))

	gwas_scores = collections.defaultdict(lambda: collections.defaultdict(lambda: 1))
	for gwas_snp in gwas_snps:
		for gwas_association in gwas_snp.evidence:
			if gwas_association.pvalue < gwas_scores[gwas_association.source][gwas_snp.snp.rsID]:
				gwas_scores[gwas_association.source][gwas_snp.snp.rsID] = gwas_association.pvalue

	functional_scores = collections.defaultdict(lambda: collections.defaultdict(int))
	snp_scores = collections.defaultdict(int)
	for gene_snp_association in association.evidence:
		for evidence in gene_snp_association.regulatory_evidence + gene_snp_association.cisregulatory_evidence:
			functional_scores[evidence.snp.rsID][evidence.source] += evidence.score
		snp_scores[gene_snp_association.snp.rsID] = gene_snp_association.score
	
	pretty_strings = []
	for ld_snp in cluster.ld_snps:
		if snp_scores[ld_snp.rsID] > 0:
			results = [ld_snp.rsID, gene_name, gene_id, ",".join(disease_names), ",".join(disease_efos), str(snp_scores[ld_snp.rsID])]
			results.append(",".join(gwas_snp.snp.rsID for gwas_snp in gwas_snps))
			for gwas_source in ['GWAS Catalog', 'GRASP', 'GWAS DB', 'Phewas Catalog']:
				results.append(",".join(str(gwas_scores[gwas_source][gwas_snp.snp.rsID]) for gwas_snp in cluster.gwas_snps))
			results += [str(functional_scores[ld_snp.rsID][functional_source]) for functional_source in ['VEP', 'Regulome', 'GTEx', 'DHS', 'Fantom5', 'PCHIC', 'Nearest']]
			pretty_strings.append("\t".join(results))

	return "\n".join(pretty_strings) + "\n"

# List of databases used
database_functions = [GWASCatalog, GWAS_DB, Phewas_Catalog, GRASP]
ld_snp_to_gene_functions = [GTEx, Fantom5, VEP, DHS, PCHIC]
snp_regulatory_functions = [Regulome] # TODO Implement and insert GERP code

if __name__ == "__main__":
	main()
