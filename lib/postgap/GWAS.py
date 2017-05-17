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
import sys

import postgap.REST
import postgap.Globals
from postgap.DataModel import *
from postgap.Utils import *

import logging
import sys

#postgap.REST.DEBUG = False


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

		logger = logging.getLogger(__name__)

		if efos is not None and len(efos) > 0:
			res = concatenate(self.query(query) for query in efos)
		else:
			res = concatenate(self.query(query) for query in diseases)

		logger.debug("\tFound %i GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in GWAS Catalog" % (len(res), ", ".join(diseases), ", ".join(efos)))

		return res

	def query(self, efo):
		logger = logging.getLogger(__name__)
		logger.info("Querying GWAS catalog for " + efo);

		server = 'http://wwwdev.ebi.ac.uk'
		url = '/gwas/beta/rest/api/efoTraits/search/findByUri?uri=%s' % (efo)
		
		#print "Querying: " + server + url;

		hash = postgap.REST.get(server, url)
		list_of_GWAS_Associations = []

		efoTraits = hash["_embedded"]["efoTraits"]

		for efoTraitHash in efoTraits:

			efoTraitLinks = efoTraitHash["_links"]
			efoTraitName  = efoTraitHash["trait"]

			logger.info("Querying Gwas rest server for SNPs associated with " + efoTraitName)

			association_rest_response = efoTraitLinks["associations"]
			association_url = association_rest_response["href"]
			association_response = postgap.REST.get(association_url, "")
			associations = association_response["_embedded"]["associations"]

			logger.info("Received " + str(len(associations)) + " associations with SNPs.")
			logger.info("Fetching SNPs and pvalues.")

			for current_association in associations:

				snp_url = current_association["_links"]["snps"]["href"]
				snp_response = postgap.REST.get(snp_url, "")
				singleNucleotidePolymorphisms = snp_response["_embedded"]["singleNucleotidePolymorphisms"]

				if (len(singleNucleotidePolymorphisms) == 0):
					sys.exit("Got no snp for a pvalue!")

				study_url = current_association["_links"]["study"]["href"]
				study_response = postgap.REST.get(study_url, "")
				pubmedId = study_response["pubmedId"]

				diseaseTrait_url = study_response["_links"]["diseaseTrait"]["href"]
				diseaseTrait_response = postgap.REST.get(diseaseTrait_url, "")
				diseaseTrait = diseaseTrait_response['trait']

				for current_snp in singleNucleotidePolymorphisms:
					logger.debug("    received association with snp rsId: " + '{:12}'.format(current_snp["rsId"]) + " with a pvalue of " + str(current_association["pvalue"]))

					list_of_GWAS_Associations.append(
						GWAS_Association(
							disease = Disease(
								name = efoTraitName,
								efo  = efo
							),
							reported_trait = diseaseTrait,
							snp     = current_snp["rsId"],
							pvalue  = current_association["pvalue"],
							source  = 'GWAS Catalog',
							study   = 'PMID' + pubmedId
						)
					)
		if len(list_of_GWAS_Associations) > 0:
			logger.info("Successfully fetched " +  str(len(list_of_GWAS_Associations)) + " SNPs and pvalues.")
		if len(list_of_GWAS_Associations) == 0:
			logger.info("Found no associated SNPs and pvalues.")
	
		return list_of_GWAS_Associations

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
		logger = logging.getLogger(__name__)

		# Temporary hack: current files have IDs, not IRIs:
		efos = [re.sub('.*/',"", efo) for efo in efos]
		
		file = open(postgap.Globals.DATABASES_DIR+"/GRASP.txt")
		res = [ self.get_association(line, diseases, efos) for line in file ]
		res = filter(lambda X: X is not None, res)

		logger.info("\tFound %i GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in GRASP" % (len(res), ", ".join(diseases), ", ".join(efos)))

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
		if items[12] in diseases or items[70] in efos:
			try:
				return GWAS_Association(
					pvalue = float(items[10]),
					snp = "rs" + items[4],
					disease = Disease(name = postgap.EFO.term(items[70]), efo = items[70]),
					reported_trait = items[12],
					source = self.display_name,
					study = items[7]
				)
			except:
				return None
		else:
			return None

class Phewas_Catalog(GWAS_source):
	display_name = "Phewas Catalog"
	logger = logging.getLogger(__name__)
	
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

		self.logger.info("\tFound %i GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in Phewas Catalog" % (len(res), ", ".join(diseases), ", ".join(efos)))

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
		# Temporary hack: current files have IDs, not IRIs:
		efos = [re.sub('.*/',"", efo) for efo in efos]
		if items[2] in diseases or items[9] in efos:
			return GWAS_Association (
				pvalue = float(items[4]),
				snp = items[1],
				disease = Disease(name = postgap.EFO.term(items[9]), efo = items[9]), 
				reported_trait = items[2],
				source = self.display_name,
				study = None
			)
		else:
			return None

class GWAS_DB(GWAS_source):
	display_name = "GWAS DB"
	logger = logging.getLogger(__name__)
	
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
		efos2 = [re.sub("_", "ID:", re.sub(".*/","", efo)) for efo in efos]
	
		res = [ self.get_association(line, diseases, efos) for line in file ]
		res = filter(lambda X: X is not None, res)

		self.logger.info("\tFound %i GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in GWAS DB" % (len(res), ", ".join(diseases), ", ".join(efos)))

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
				disease = Disease(name = postgap.EFO.term(items[6]), efo = items[6]),
				reported_trait = items[5].decode('latin1'),
				source = self.display_name,
				study = items[4]
			)
		else:
			return None

sources = GWAS_source.__subclasses__()
