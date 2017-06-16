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
import cPickle as pickle
import collections

from postgap.DataModel import *
from postgap.Utils import *
import REST
import postgap.Globals 
import BedTools
import Ensembl_lookup
import requests
import logging

VEP_impact_to_score = {
	'HIGH': 4,
	'MEDIUM': 3,
	'LOW': 2,
	'MODIFIER': 1,
	'MODERATE': 1
}

class Cisreg_source(object):
	logger = logging.getLogger(__name__)
	def run(self, snps, tissues):
		assert False, "This stub should be defined"

class GTEx(Cisreg_source):
	display_name = "GTEx"

	
	def run(self, snps, tissues):
		"""

			Returns all genes associated to a set of SNPs in GTEx
			Args:
			* [ SNP ]
			* [ string ] (tissues)
			Returntype: [ Cisregulatory_Evidence ]

		"""
		# Find all genes with 1Mb
		start = min(snp.pos for snp in snps)
		end = max(snp.pos for snp in snps)
		chrom = snps[0].chrom
		

		server = 'http://grch37.rest.ensembl.org'
		ext = '/overlap/region/%s/%s:%i-%i?feature=gene;content-type=application/json' % (postgap.Globals.SPECIES, chrom, max(0, start - 1e6), end + 1e6)
		genes = [ Gene(
				name = gene['external_name'],
				id = gene['id'],
				chrom = gene['seq_region_name'],
				tss = int(gene['start']) if gene['strand'] > 0 else int(gene['end']),
				biotype = gene['biotype']
			)
			for gene in postgap.REST.get(server, ext)
		]
		if len(genes) < len(snps):
			snp_hash = dict( (snp.rsID, snp) for snp in snps)
			res = concatenate((self.gene(gene, tissues, snp_hash) for gene in genes))
		else:
			res = concatenate((self.snp(snp, tissues) for snp in snps))

		logger = logging.getLogger(__name__)
		logger.info("\tFound %i interactions in GTEx" % (len(res)))
		print "\tFound %i interactions in GTEx" % (len(res))
		#raise Exception("Being called!")
		return res

	def gene(self, gene, tissues, snp_hash):
		"""

			Returns all SNPs associated to a gene in GTEx
			Args:
			* Gene
			* [ string ] (tissues)
			* { rsID: SNP }
			Returntype: [ Cisregulatory_Evidence ]

		"""
		res = concatenate(self.gene_tissue(gene, tissue, snp_hash) for tissue in tissues)

		self.logger.info("\tFound %i SNPs associated to gene %s in GTEx" % (len(res), gene.id))

		return res

	def gene_tissue(self, gene, tissue, snp_hash):
		"""

			Returns all SNPs associated to a gene in GTEx in a given tissue
			Args:
			* Gene
			* string, tissue
			* { rsID: SNP }
			Returntype: [ Cisregulatory_Evidence ]

		"""

		server = "http://rest.ensembl.org"
		ext = "/eqtl/id/%s/%s?content-type=application/json;statistic=p-value;tissue=%s" % ('homo_sapiens', gene.id, tissue);
		try:

			eQTLs = postgap.REST.get(server, ext)

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
					source = self.display_name,
					study = None,
					info = None,
					pvalue = eQTL['value'],
					beta = None
				)
				for eQTL in eQTLs 
				if eQTL['snp'] in snp_hash
				if eQTL['value'] < 2.5e-5 
			]

			self.logger.info("\tFound %i SNPs associated to gene %s in tissue %s in GTEx" % (len(res), gene.id, tissue))

			return res
		except Exception as e:
			self.logger.warning("Got exception when querying %s%s" % (server, ext))
			self.logger.warning("The exception is %s" % (e))
			self.logger.warning("Returning 'None' and pretending this didn't happen.")
			return None

	def gene_tissue_betas(self, gene, tissue, snp_hash):
		"""

			Returns all SNPs associated to a gene in GTEx in a given tissue
			Args:
			* Gene
			* string, tissue
			* { rsID: SNP }
			Returntype: [ Cisregulatory_Evidence ]

		"""
		server = "http://rest.ensembl.org"
		ext = "/eqtl/id/%s/%s?content-type=application/json;statistic=beta;tissue=%s" % ('homo_sapiens', gene.id, tissue);
		try:
			eQTLs = postgap.REST.get(server, ext)

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
					source = self.display_name,
					study = None,
					info = None,
					pvalue = None,
					beta = eQTL['value']
				)
				for eQTL in eQTLs 
				if eQTL['snp'] in snp_hash
			]

			self.logger.info("\tFound %i SNPs with betas associated to gene %s in tissue %s in GTEx" % (len(res), gene.id, tissue))

			return res
		except Exception as e:
			self.logger.warning("Got exception when querying %s%s" % (server, ext))
			self.logger.warning("The exception is %s" % (e))
			self.logger.warning("Returning 'None' and pretending this didn't happen.")
			return None

	def snp(self, snp, tissues):
		"""

			Returns all genes associated to a snp in GTEx
			Args:
			* SNP
			* [ string ] (tissues)
			Returntype: [ Cisregulatory_Evidence ]

		"""
		res = concatenate(self.snp_tissue(snp, tissue) for tissue in tissues)

		self.logger.info("\tFound %i genes associated to snp %s in GTEx" % (len(res), snp.rsID))

		return res

	def snp_tissue(self, snp, tissue):
		"""

			Returns all genes associated to a SNP in GTEx in a given tissue
			Args:
			* SNP
			* string, tissues
			Returntype: [ Cisregulatory_Evidence ]

		"""


		server = "http://rest.ensembl.org"
		ext = "/eqtl/variant_name/%s/%s?content-type=application/json;statistic=p-value;tissue=%s" % ('homo_sapiens', snp.rsID, tissue);
		try:
			eQTLs = postgap.REST.get(server, ext)

			'''
				Example return object:
				[
					{
						minus_log10_p_value: 1.47569690641653,
						value: 0.0334428355738418,
						gene: "ENSG00000162627"
					},
				]
			'''

			res = [
				Cisregulatory_Evidence(
					snp = snp,
					gene = postgap.Ensembl_lookup.get_ensembl_gene(eQTL['gene']),
					tissue = tissue,
					score = 1,
					source = self.display_name,
					study = None,
					info = None
				)
				for eQTL in eQTLs
				if eQTL['value'] < 2.5e-5 
			]

			self.logger.info("\tFound %i genes associated the SNP %s in tissue %s in GTEx" % (len(res), snp.rsID, tissue))

			return res
		except:
			return None

class VEP(Cisreg_source):
	display_name = "VEP"

	
	def run(self, snps, tissues):
		"""

			Returns all genes associated to a set of SNPs in VEP
			Args:
			* [ SNP ]
			* [ string ] (tissues)
			Returntype: [ Regulatory_Evidence ]

		"""
		list = concatenate(self.get(chunk) for chunk in chunks(snps, 199))
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

		snp_hash = dict( (snp.rsID, snp) for snp in snps)
		transcript_consequences = filter(lambda X: 'transcript_consequences' in X, list)
		res = []
		for hit in transcript_consequences:
			for consequence in hit['transcript_consequences']:
				res.append(Cisregulatory_Evidence(
					snp = snp_hash[hit['input']],
					gene = postgap.Ensembl_lookup.get_ensembl_gene(consequence['gene_id']),
					score = VEP_impact_to_score[consequence['impact']],
					source = self.display_name,
					study = None,
					tissue = None,
					info = {
						'consequence_terms': consequence['consequence_terms'], 
					}
				))

		self.logger.info("\tFound %i interactions in VEP" % (len(res)))

		return res

	def remove_none_elements(self, list):
		return filter(self.exists, list)

	def exists(self, it):
		return (it is not None)

	def get(self, chunk_param):
		"""

			Queries Ensembl servers. Recursively breaks down query if error code 400 
			("Bad request") is returned
			Args:
			* [ SNP ]
			Returntype: [ Regulatory_Evidence ]

		"""

		# 
		chunk = self.remove_none_elements(chunk_param)

		if len(chunk) == 0:
			return []

		try:
			server = "http://grch37.rest.ensembl.org"
			ext = "/vep/%s/id" % (postgap.Globals.SPECIES)
			return postgap.REST.get(server, ext, data = {"ids" : [snp.rsID for snp in chunk]})

		except requests.exceptions.HTTPError as error:
			#if error.response.status_code == 400 or error.response.status_code == 504:
			if \
				error.response.status_code == requests.codes.bad_request \
				or error.response.status_code == requests.codes.gateway_timeout:

				if len(chunk) == 1:
					return []
				else:
					return self.get(chunk[:len(chunk)/2]) + self.get(chunk[len(chunk)/2:])
			raise
				

class Fantom5(Cisreg_source):
	display_name = "Fantom5"

	
	def run(self, snps, tissues):
		"""

			Returns all genes associated to a set of SNPs in Fantom5
			Args:
			* [ SNP ]
			* [ string ] (tissues)
			Returntype: [ Regulatory_Evidence ]

		"""
		
		self.logger.info("\tSearching for overlaps from %i SNPs to Fantom5" % len(snps))
		
		intersection = postgap.BedTools.overlap_snps_to_bed(snps, postgap.Globals.DATABASES_DIR + "/Fantom5.bed")
		fdr_model = pickle.load(open(postgap.Globals.DATABASES_DIR + "/Fantom5.fdrs"))
		snp_hash = dict( (snp.rsID, snp) for snp in snps)
		hits  = filter(lambda X: X is not None, map(lambda X: self.get_evidence(X, fdr_model, snp_hash), intersection))

		self.logger.info("\tFound %i overlaps with %i SNPs to Fantom5" % (len(hits), len(snps)))

		res = filter(lambda X: X.score, hits)

		self.logger.info("\tFound %i interactions in Fantom5" % (len(res)))

		return res

	def get_evidence(self, feature, fdr_model, snp_hash):
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
		gene = postgap.Ensembl_lookup.get_gene(feature[3])
		if gene is None:
			return None
		snp = snp_hash[feature[8]]
		score = STOPGAP_FDR(snp, gene, fdr_model)

		return Cisregulatory_Evidence(
				snp = snp,
				gene = gene,
				source = self.display_name,
				score = score,
				study = None,
				tissue = None,
				info = None
			)

class DHS(Cisreg_source):
	display_name = "DHS"

	
	def run (self, snps, tissues):
		"""

			Returns all genes associated to a set of SNPs in DHS
			Args:
			* [ SNP ]
			* [ string ] (tissues)
			Returntype: [ Regulatory_Evidence ]

		"""
		self.logger.info("\tSearching for gene associations in DHS")
		
		intersection = postgap.BedTools.overlap_snps_to_bed(snps, postgap.Globals.DATABASES_DIR + "/DHS.bed")
		fdr_model = pickle.load(open(postgap.Globals.DATABASES_DIR+"/DHS.fdrs"))
		snp_hash = dict( (snp.rsID, snp) for snp in snps)
		res = filter (lambda X: X is not None and X.score, (self.get_evidence(feature, fdr_model, snp_hash) for feature in intersection))

		self.logger.info("\tFound %i gene associations in DHS" % len(res))

		return res

	def get_evidence(self, feature, fdr_model, snp_hash):
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
		gene = postgap.Ensembl_lookup.get_gene(feature[3])
		if gene is None:
			return None
		snp = snp_hash[feature[8]]
		score = STOPGAP_FDR(snp, gene, fdr_model)

		return Cisregulatory_Evidence(
			snp = snp,
			gene = gene,
			score = score,
			source = self.display_name,
			study = None,
			tissue = None,
			info = None
		)

def get_fdr_model(filename):
	return pickle.load(open(filename))

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

	if distance > fdr_model['MAX_DISTANCE']:
		return 0

	bin = int(distance / fdr_model['BIN_WIDTH'])

	if bin not in fdr_model['FDR']:
		return 0

	FDR = fdr_model['FDR'][bin]

	if FDR is None:
		return 0
	elif FDR < .6:
		return 2
	elif FDR < .85:
		return 1
	else:
		return 0

class PCHIC(Cisreg_source):
	display_name = "PCHiC"

	
	def run(self, snps, tissues):
		"""

			Returns all genes associated to a set of SNPs in PCHIC
			Args:
			* [ SNP ]
			* [ string ] (tissues)
			Returntype: [ Regulatory_Evidence ]

		"""
		
		self.logger.info("\tSearching for gene associations in PCHIC")
		
		intersection = postgap.BedTools.overlap_snps_to_bed(snps, postgap.Globals.DATABASES_DIR + "/pchic.bed")
		snp_hash = dict( (snp.rsID, snp) for snp in snps)
		res = filter (lambda X: X is not None and X.score, (self.get_evidence(feature, snp_hash) for feature in intersection))

		self.logger.info("\tFound %i gene associations in PCHIC" % len(res))

		return res

	def get_evidence(self, feature, snp_hash):
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
		gene = postgap.Ensembl_lookup.get_gene(feature[4])
		if gene is None:
			return None
		snp = snp_hash[feature[8]]

		return Cisregulatory_Evidence(
			snp = snp,
			gene = gene,
			score = 1,
			source = self.display_name,
			study = None,
			tissue = None,
			info = None
		)

class nearest_gene(Cisreg_source):
	display_name = 'Nearest'
	
	def run (self, snps, tissues):
		"""

			Return nearest gene to SNP
			Args:
			* [ SNP ]
			Returntype: dict(rsID => Gene)

		"""
		snps = list(snps)
		bed = postgap.Globals.DATABASES_DIR + "/Ensembl_TSSs.bed"
		res = postgap.BedTools.closest(snps, bed)
		snp_hash = dict((snp.rsID, snp) for snp in snps)

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
		return [ Cisregulatory_Evidence(
				gene = postgap.Ensembl_lookup.get_gene(row[7]),
				snp = snp_hash[row[3]],
				score = 1,
				source = self.display_name,
				study = None,
				tissue = None,
				info = None
				)
			for row in res ]

sources = Cisreg_source.__subclasses__()
