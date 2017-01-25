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
import re
import requests
import json

import postgap.REST
import postgap.Globals
from postgap.DataModel import *
from postgap.Utils import *

class GWAS_source(object):
	def run(self, diseases, efos):
		"""

			Returns all GWAS SNPs associated to a disease in this source
			Args:
			* [ string ] (trait descriptions)
			* [ string ] (trait EFO identifiers)
			Returntype: [ GWAS_Association ]

		"""
		assert False, "This stub should be defined"

class GWASCatalog(GWAS_source):
	display_name = 'GWAS Catalog'

	def run(self, diseases, efos):
		"""

			Returns all GWAS SNPs associated to a disease in GWAS Catalog
			Args:
			* [ string ] (trait descriptions)
			* [ string ] (trait EFO identifiers)
			Returntype: [ GWAS_Association ]

		"""
		if efos is not None and len(efos) > 0:
			res = concatenate(self.query(query) for query in efos)
		else:
			res = concatenate(self.query(query) for query in diseases)

		if postgap.Globals.DEBUG:
			print "\tFound %i GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in GWAS Catalog" % (len(res), ", ".join(diseases), ", ".join(efos))

		return res

	def query(self, term):
		server = 'http://www.ebi.ac.uk'
		url_term = re.sub(" ", "%20", term)
		ext1 = '/gwas/api/search/moreresults?q="%s"&max=0&facet=association&pvalfilter=&orfilter=&betafilter=&datefilter=&sort=' % (url_term)
		hash = postgap.REST.get(server, ext1)
		try:
			count = int(hash['response']['numFound'])
		except:
			print "Failed on %s%s" % (server, ext1)
		ext2 = '/gwas/api/search/moreresults?q="%s"&max=%i&facet=association&pvalfilter=&orfilter=&betafilter=&datefilter=&sort=' % (url_term, count)
		hash = postgap.REST.get(server, ext2)
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
					source = self.display_name,
					study = "PMID" + hit['pubmedId']
				)
				for hit in hits
			]
		except:
			print "Failed on %s%s" % (server, ext2)

class GRASP(GWAS_source):
	display_name = "GRASP"
	def run(self, diseases, efos):
		"""

			Returns all GWAS SNPs associated to a disease in GRASP
			Args:
			* [ string ] (trait descriptions)
			* [ string ] (trait EFO identifiers)
			Returntype: [ GWAS_Association ]

		"""
		file = open(postgap.Globals.DATABASES_DIR+"/GRASP.txt")
		res = [ self.get_association(line, diseases, efos) for line in file ]
		res = filter(lambda X: X is not None, res)

		if postgap.Globals.DEBUG:
			print "\tFound %i GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in GRASP" % (len(res), ", ".join(diseases), ", ".join(efos))

		return res

	def get_association(self, line, diseases, efos):
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
		if items[11] in diseases or items[70] in efos:
			return GWAS_Association(
				pvalue = float(items[10]),
				snp = "rs" + items[4],
				disease = Disease(name = items[11], efo = items[70]),
				source = self.display_name,
				study = items[7]
			)
		else:
			return None

class Phewas_Catalog(GWAS_source):
	display_name = "Phewas Catalog"
	def run(self, diseases, efos):
		"""

			Returns all GWAS SNPs associated to a disease in PhewasCatalog
			Args:
			* [ string ] (trait descriptions)
			* [ string ] (trait EFO identifiers)
			Returntype: [ GWAS_Association ]

		"""
		file = open(postgap.Globals.DATABASES_DIR+"/Phewas_Catalog.txt")
		res = [ self.get_association(line, diseases, efos) for line in file ]
		res = filter(lambda X: X is not None, res)

		if postgap.Globals.DEBUG:
			print "\tFound %i GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in Phewas Catalog" % (len(res), ", ".join(diseases), ", ".join(efos))

		return res

	def get_association(self, line, diseases, efos):
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
				source = self.display_name,
				study = None
			)
		else:
			return None

class GWAS_DB(GWAS_source):
	display_name = "GWAS DB"
	def run(self, diseases, efos):
		"""

			Returns all GWAS SNPs associated to a disease in GWAS_DB
			Args:
			* [ string ] (trait descriptions)
			* [ string ] (trait EFO identifiers)
			Returntype: [ GWAS_Association ]

		"""
		file = open(postgap.Globals.DATABASES_DIR+"/GWAS_DB.txt")

		# Note the format of the EFO strings is modified in this file format, so we need to change the queries
		efos2 = [re.sub("_", "ID:", efo) for efo in efos]
		res = [ self.get_association(line, diseases, efos) for line in file ]
		res = filter(lambda X: X is not None, res)

		if postgap.Globals.DEBUG:
			print "\tFound %i GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in GWAS DB" % (len(res), ", ".join(diseases), ", ".join(efos))

		return res

	def get_association(self, line, diseases, efos):
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
				source = self.display_name,
				study = items[4]
			)
		else:
			return None

sources = GWAS_source.__subclasses__()
