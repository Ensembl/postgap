#! /usr/bin/env python

"""

Copyright [1999-2018] EMBL-European Bioinformatics Institute

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
from postgap.REST import Variation400error
import postgap.REST
import postgap.Globals 
import BedTools
import Ensembl_lookup
import requests
import logging
import sqlite3
import h5py
import sys
import numpy
import postgap.FinemapIntegration


VEP_impact_to_score = {
	'HIGH': 4,
	'MEDIUM': 3,
	'LOW': 2,
	'MODIFIER': 1,
	'MODERATE': 1
}

class Cisreg_source(object):
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

		res = concatenate(map(self.snp, snps))
		logging.info("\tFound %i interactions in GTEx" % (len(res)))
		return res

	def snp(self, snp):
		"""

			Returns all SNPs associated to a SNP in GTEx
			Args:
			* string (rsID)
			Returntype: [ Cisregulatory_Evidence ]

		"""
		if postgap.Globals.GTEx_path is None:
			res = self._snp_rest(snp)
		else:
			res = self._snp_hdf5(snp)

		number_of_associated_genes = 0
		if res is not None:
			number_of_associated_genes = len(res)
			

		logging.info(
			"\tFound %i genes associated the SNP %s in GTEx" % (number_of_associated_genes, snp.rsID))

		return res


	def _snp_rest(self, snp):
		"""

			Returns all SNPs associated to a SNP in GTEx
			Args:
			* string (rsID)
			Returntype: [ Cisregulatory_Evidence ]

		"""
		cisreg_with_pvalues = self._snp_pvalues(snp)
		
		if cisreg_with_pvalues is None:
			return []


		if not postgap.Globals.PERFORM_BAYESIAN:
			return cisreg_with_pvalues


		if len(cisreg_with_pvalues) == 0:
			# Empty list
			return cisreg_with_pvalues
		
		cisreg_betas = self._snp_betas(snp)
		
		# Where there are pvalues, there must be betas
		if cisreg_betas is None:
			logging.warning("Got exception in _snp_rest")
			logging.warning("The exception is cisreg_betas is None")
			logging.warning("Returning 'None' and pretending this didn't happen.")
			return None
			#raise Exception
		
		if len(cisreg_betas) == 0:
			logging.warning("Got exception in _snp_rest")
			logging.warning("The exception is the length of cisreg_betas is 0")
			logging.warning("Returning 'None' and pretending this didn't happen.")
			return None
			#raise Exception
		
		# Match them up:
		
		combined_cisreg_evidence_list = []
		
		for cisreg_with_pvalue in cisreg_with_pvalues:
			try:
				matching_cisreg_betas = filter(lambda X: X.gene.id == cisreg_with_pvalue.gene.id and X.tissue == cisreg_with_pvalue.tissue, cisreg_betas)

			except Exception as e:
				logging.warning("Got exception matching_cisreg_betas in _snp_rest")
				logging.warning("The exception is %s" % (e))
				continue
			
			if len(matching_cisreg_betas) != 1:
				logging.warning("Got exception in _snp_rest")
				logging.warning("The exception is matching_cisreg_betas != 1; (matching_cisreg_betas = %i)" % len(matching_cisreg_betas))
				continue
				#raise Exception
			
			cisreg_with_beta = matching_cisreg_betas[0]

			if postgap.Globals.PERFORM_BAYESIAN:
				z_score = postgap.FinemapIntegration.z_score_from_pvalue(cisreg_with_pvalue.pvalue,cisreg_with_beta.beta),
			else:
				z_score = None

			assert cisreg_with_beta.beta is not None

			combined_cisreg_evidence = Cisregulatory_Evidence(
				snp    = cisreg_with_pvalue.snp,
				gene   = cisreg_with_pvalue.gene,
				tissue = cisreg_with_pvalue.tissue,
				score  = 1 - cisreg_with_pvalue.pvalue,
				source = cisreg_with_pvalue.source,
				study  = 'GTEx',
				info   = None,
				z_score = z_score, 
				pvalue = cisreg_with_pvalue.pvalue,
				beta   = cisreg_with_beta.beta
			)
			combined_cisreg_evidence_list.append(combined_cisreg_evidence)
		
		return combined_cisreg_evidence_list

	def _snp_betas(self, snp):
		"""

			Returns all genes associated to a snp in GTEx
			Args:
			* SNP
			Returntype: [ Cisregulatory_Evidence ]

		"""
		try:
			server = "http://rest.ensembl.org"
			ext = "/eqtl/variant_name/%s/%s?content-type=application/json;statistic=beta" % ('homo_sapiens', snp.rsID);

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

			return [
				Cisregulatory_Evidence(
					snp = snp,
					gene = postgap.Ensembl_lookup.get_ensembl_gene(eQTL['gene']),
					tissue = eQTL['tissue'],
					score = 1,
					source = self.display_name,
					study = None,
					info = None,
					z_score = None,
					pvalue = None,
					beta = float(eQTL['value'])
				)
				for eQTL in eQTLs
			]

		except Exception as e:
			logging.warning("Got exception when quering _snp_betas")
			logging.warning("The exception is %s" % (e))
			logging.warning("Returning 'None' and pretending this didn't happen.")
			return None

	def _snp_pvalues(self, snp):
		"""

			Returns all genes associated to a snp in GTEx
			Args:
			* SNP
			Returntype: [ Cisregulatory_Evidence ]

		"""
		try:
			server = "http://rest.ensembl.org"
			ext = "/eqtl/variant_name/%s/%s?content-type=application/json;statistic=p-value" % (
				'homo_sapiens', snp.rsID);

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
			return [
				Cisregulatory_Evidence(
					snp = snp,
					gene = postgap.Ensembl_lookup.get_ensembl_gene(eQTL['gene']),
					tissue = eQTL['tissue'],
					score = 1 - float(eQTL['value']),
					source = self.display_name,
					study = None,
					info = None,
					pvalue = float(eQTL['value']),
					beta = None,
					z_score = None
				)
				for eQTL in eQTLs
			]

		except Exception, e:
			logging.warning("Got exception when quering _snp_pvalues")
			logging.warning("The exception is %s" % (e))
			logging.warning("Returning 'None' and pretending this didn't happen.")
			return None

	def _snp_hdf5(self, snp):
		"""

			Returns all SNPs associated to a SNP in GTEx
			Args:
			* SNP
			Returntype: [ Cisregulatory_Evidence ]

		"""
		try:

			snp_cursor = get_sqlite_values(['hdf5_index', 'external_id'], 'snp', 'external_id', [snp.rsID])
			if snp_cursor is None:
				raise ValueError('SNP: empty cursor')

			snp_val = snp_cursor.fetchone()

			if snp_val is None:
				raise ValueError('SNP: empty fetch')

			hdf5_snp_index = snp_val[0]
			hdf5_snp_index = hdf5_snp_index  - 1


			with h5py.File(postgap.Globals.GTEx_path, 'r') as hdf5_file:
				tissue_array_names = {k: (''.join(chr(i) for i in hdf5_file.get('dim_labels/1')[k]).rstrip('\0'))
									  for k in range(0, len(hdf5_file.get('dim_labels/1')))}

				hdf5_gene_boundaries = hdf5_file.get('boundaries/3')[hdf5_snp_index]


				beta = numpy.array([hdf5_file['matrix'][0, :, int(hdf5_gene_boundaries[0][0]):int(hdf5_gene_boundaries[0][1]) + 1, hdf5_snp_index]])
				p_val = numpy.array([hdf5_file['matrix'][1, :, int(hdf5_gene_boundaries[0][0]):int(hdf5_gene_boundaries[0][1]) + 1, hdf5_snp_index]])

				p_val_index = numpy.where(p_val > 0)

				gene_range = range(int(hdf5_gene_boundaries[0][0]), int(hdf5_gene_boundaries[0][1]) + 1)
				gene_range_filtered = [gene_range[k] for k in p_val_index[2]]
				gene_index_list = set(p_val_index[2])
				gene_list_filtered = [gene_range[k] for k in gene_index_list]
				gene_array_names = {k:(''.join(chr(i) for i in hdf5_file.get('dim_labels/2')[k])).rstrip('\0') for k in gene_list_filtered}


			res = []
			for k in range(len(p_val_index[0])):
				pvalue = p_val[p_val_index[0][k]][p_val_index[1][k]][p_val_index[2][k]]
				
				if pvalue > 1 or pvalue < 0:
					continue

				try:
					beta_val = beta[p_val_index[0][k]][p_val_index[1][k]][p_val_index[2][k]]
				except Exception as e:
					logging.warning("Got exception when quering getting BETA value from snp_hdf5")
					logging.warning("The exception is %s" % (e))					
					continue

				if postgap.Globals.PERFORM_BAYESIAN:
					z_score = postgap.FinemapIntegration.z_score_from_pvalue(pvalue,beta_val),
				else:
					z_score = None


				res.append([
					Cisregulatory_Evidence(
						snp = snp,
						gene= postgap.Ensembl_lookup.get_ensembl_gene(gene_array_names[gene_range_filtered[k]]),
						tissue = tissue_array_names[p_val_index[1][k]],
						score = 1 - float(pvalue),
						source = self.display_name,
						study = None,
						info = None,
						z_score = z_score,
						pvalue = pvalue,
						beta = float(beta_val)
					)
				])

			return res

		except Exception as e:
			logging.warning("Got exception when quering snp_hdf5")
			logging.warning("The exception is %s" % (e))
			logging.warning("Returning 'None' and pretending this didn't happen.")
			return []

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

		list = concatenate(self.get(chunk) for chunk in chunks(snps, 1))
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
					},
					z_score = None,
					pvalue = None,
					beta = None
				))

		logging.info("\tFound %i interactions in VEP" % (len(res)))

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
		except Variation400error as error:
			if len(chunk) == 1:
				return []
			else:
				return self.get(chunk[:len(chunk)/2]) + self.get(chunk[len(chunk)/2:])
				

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

		logging.info("\tSearching for overlaps from %i SNPs to Fantom5" % len(snps))
		
		intersection = postgap.BedTools.overlap_snps_to_bed(snps, postgap.Globals.DATABASES_DIR + "/Fantom5.bed")
		fdr_model = pickle.load(open(postgap.Globals.DATABASES_DIR + "/Fantom5.fdrs"))
		snp_hash = dict( (snp.rsID, snp) for snp in snps)
		hits  = filter(lambda X: X is not None, map(lambda X: self.get_evidence(X, fdr_model, snp_hash), intersection))

		logging.info("\tFound %i overlaps with %i SNPs to Fantom5" % (len(hits), len(snps)))

		res = filter(lambda X: X.score, hits)

		logging.info("\tFound %i interactions in Fantom5" % (len(res)))

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
				info = None,
				z_score = None,
				pvalue = None,
				beta = None
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

		logging.info("\tSearching for gene associations in DHS")
		
		intersection = postgap.BedTools.overlap_snps_to_bed(snps, postgap.Globals.DATABASES_DIR + "/DHS.bed")
		fdr_model = pickle.load(open(postgap.Globals.DATABASES_DIR+"/DHS.fdrs"))
		snp_hash = dict( (snp.rsID, snp) for snp in snps)
		res = filter (lambda X: X is not None and X.score, (self.get_evidence(feature, fdr_model, snp_hash) for feature in intersection))

		logging.info("\tFound %i gene associations in DHS" % len(res))

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
			info = None,
			z_score = None,
			pvalue = None,
			beta = None
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

	return 1 - FDR

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

		logging.info("\tSearching for gene associations in PCHIC")
		
		intersection = postgap.BedTools.overlap_snps_to_bed(snps, postgap.Globals.DATABASES_DIR + "/pchic.bed")
		snp_hash = dict( (snp.rsID, snp) for snp in snps)
		res = filter (lambda X: X is not None and X.score, (self.get_evidence(feature, snp_hash) for feature in intersection))

		logging.info("\tFound %i gene associations in PCHIC" % len(res))

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
			score = float(feature[3]) / 645.0974, # Divide by max value in dataset
			source = self.display_name,
			study = None,
			tissue = None,
			info = None,
			z_score = None,
			pvalue = None,
			beta = None
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
				info = None,
				z_score = None,
				pvalue = None,
				beta = None
				)
			for row in res ]


def get_filtered_subclasses(subclasses_filter):
	return [subclass for subclass in sources if subclass.display_name in subclasses_filter]


def get_sqlite_values(columns, dimension, filtername, args):
	try:
		sql = "SELECT {fld}".format(fld=', '.join(columns)) + " FROM %s" % dimension
		if filtername is not None:
			sql += " WHERE %s" % filtername + " in ({seq})".format(seq=', '.join(map("'{0}'".format, args)))

		db = sqlite3.connect(postgap.Globals.SQLite_connection)

		cursor = db.cursor().execute(sql)
		db.close

	except Exception as e:
		logging.warning("Got exception in get_sqlite_values")
		logging.warning("The exception is %s" % (e))
		logging.warning("Returning 'None' and pretending this didn't happen.")
		return None

	return cursor

sources = Cisreg_source.__subclasses__()
