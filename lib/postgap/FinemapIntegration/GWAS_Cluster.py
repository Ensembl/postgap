#!/usr/bin/env python

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
import sys
import json
import postgap.LD
import postgap.Globals
from postgap.DataModel import GWAS_Cluster
import logging

import collections

class ZScoreComputationException(Exception):
	pass

class FinemapFailedException(Exception):
	pass

class IntegrityCheckException(Exception):
	pass

def compute_finemap_posteriors(lead_snps, ld_snps, cluster_name="Unnamed cluster"):
	"""
		Computes finemap posteriors for a gwas cluster.
		
		Returns the finemap_posteriors returned by the as computed by the 
		stochastic search in finemap.
	"""
	
	(snps_with_z_scores, r2_array, SNP_ids) = compute_approximated_z_scores(
		lead_snp_like_objects = lead_snps,
		ld_snps               = ld_snps
	)
	
	from postgap.Summarisers import summarise
	logging.info("The snps with zscores, that will be used by finemap, are:\n" + summarise(snps_with_z_scores))
	
	z_score_list = []
	for snp_with_z_score in snps_with_z_scores:
		
		z_score = None
		
		from postgap.DataModel import GWAS_SNP, SNP, Cisregulatory_Evidence
		
		# If this is a GWAS_SNP
		if type(snp_with_z_score) is GWAS_SNP:

			# and it has a z score, then use that
			if snp_with_z_score.z_score is not None:
				z_score = snp_with_z_score.z_score
			
			# Otherwise use the approximated z score from ld.
			else:
				z_score = snp_with_z_score.snp.approximated_zscore
		
		# If it is only a snp, then the approximated z score is used.
		if type(snp_with_z_score) is SNP:
			z_score = snp_with_z_score.approximated_zscore
		
		if type(snp_with_z_score) is Cisregulatory_Evidence:
			z_score = snp_with_z_score.z_score
		
		# z_score should have some value by now
		if z_score is None:
			raise Exception("No z score for:\n" + summarise(snp_with_z_score))
		
		z_score_list.append(z_score)
	

	logging.info("Running finemap")
	
	#BUG: approximated_gwas_zscores should only be used, if the pvalues are unknown.
	
	kstart = 1
	kmax   = 1
	max_iter = "Not used when kstart == kmax"

	from postgap.Finemap import finemap
	finemap_posteriors = finemap(
		labels       = SNP_ids,
		z_scores     = z_score_list,
		cov_matrix   = r2_array,
		n            = len(r2_array),
		kstart       = kstart,
		kmax         = kmax,
		max_iter     = max_iter,
		prior        = "independence",
		sample_label = cluster_name
	)
	logging.info("Done running finemap")
	
	if finemap_posteriors is None:
		raise Exception ("ERROR: Didn't get posteriors!")

	return finemap_posteriors

def compute_z_scores_for_snp_like_object(snp_like_objects):
	
	logger = logging.getLogger(__name__)
	
	from postgap.FinemapIntegration.GWAS_SNP import compute_z_score_for_snp_like_type
	
	snp_like_object_with_z_scores = []
	
	for snp_like_object in snp_like_objects:
		
		# This may be "None"
		z_score = compute_z_score_for_snp_like_type(snp_like_object)
		
		from postgap.DataModel import GWAS_SNP, SNP, Cisregulatory_Evidence
		
		if type(snp_like_object) is GWAS_SNP:
		
			gwas_snp_with_z_score = GWAS_SNP(
				snp       = snp_like_object.snp,
				pvalue    = snp_like_object.pvalue,
				z_score   = z_score,
				evidence  = snp_like_object.evidence,
			)
			snp_like_object_with_z_scores.append(gwas_snp_with_z_score)
			continue
		
		if type(snp_like_object) is Cisregulatory_Evidence:
			
			cisregulatory_evidence_with_zscore = Cisregulatory_Evidence(
				snp     = snp_like_object.snp,
				gene    = snp_like_object.gene,
				score   = snp_like_object.score,
				source  = snp_like_object.source,
				study   = snp_like_object.study,
				tissue  = snp_like_object.tissue,
				info    = snp_like_object.info,
				z_score = z_score,
				pvalue  = snp_like_object.pvalue,
				beta    = snp_like_object.beta, 
			)
			snp_like_object_with_z_scores.append(cisregulatory_evidence_with_zscore)
			continue
		
		raise Exception("Not implemented for type " + type(snp_like_object))
	
	return snp_like_object_with_z_scores

def compute_approximated_z_scores(lead_snp_like_objects, ld_snps):
	
	from postgap.DataModel import GWAS_SNP, SNP, Cisregulatory_Evidence
	
	all_ld_snps_as_snp_tuples = []
	
	# Some of the ld_snps may be of type GWAS_SNP as well. This can happen, if
	# they were included in the analysis because of their ld value, but there
	# was gwas data available for them and so they have a pvalue and are in 
	# this as GWAS_SNPs.
	#
	#for snp in lead_snp_like_objects + ld_snps:
	for snp in ld_snps:
		
		if type(snp) is GWAS_SNP:
			all_ld_snps_as_snp_tuples.append(snp.snp)
			continue
		
		if type(snp) is SNP:
			all_ld_snps_as_snp_tuples.append(snp)
			continue

		if type(snp) is Cisregulatory_Evidence:
			all_ld_snps_as_snp_tuples.append(snp.snp)
			continue
		
		raise Exception("Unknown type: " + str(type(snp)))
	
	logger = logging.getLogger(__name__)
	
	# LD measures the deviation from the expectation of non-association.
	(SNP_ids, r2_array) = postgap.LD.get_pairwise_ld(all_ld_snps_as_snp_tuples)

	ld_cluster_size = len(ld_snps)

	if ld_cluster_size==0:
		raise Exception("The cluster size was zero!")

	logger.info("SNPs:")
	logger.info("\n" + stringify_vector(SNP_ids))

	logger.info("Matrix of pairwise linkage disequilibria:")
	logger.info("\n" + stringify_matrix(r2_array))
	
	# Only the lead snps will be used to compute the approximated gwas scores.
	# The gwas snps that were included from files, because of their ld are not
	# used to compute the approximated zscores.
	#
	z_scores_for_lead_snp_like_objects = compute_z_scores_for_snp_like_object(lead_snp_like_objects)
	
	# Collect all lead gwas snps for which z score could be computed
	lead_snp_like_objects_with_z_scores = []
	
	for snp_like_object in z_scores_for_lead_snp_like_objects:
		
		# All snp_like objects that could be lead snps are
		#
		# - GWAS_SNP and
		# - Cisregulatory_Evidence nd
		#
		# have a z_score attribute:
		#
		if not(snp_like_object.z_score is None):
			lead_snp_like_objects_with_z_scores.append(snp_like_object)
	
	if len(lead_snp_like_objects_with_z_scores) == 0:
		raise ZScoreComputationException("Couldn't compute a z score for any of the lead snps.")

	
	# This loop is just a healthcheck.
	#
	# Every lead snp should be in the SNP_id list.
	#
	for lead_snp_like_object_with_z_score in lead_snp_like_objects_with_z_scores:
		
		found_in_list = None
		try:
			SNP_ids.index(lead_snp_like_object_with_z_score.snp.rsID)
			found_in_list = True
		except ValueError:
			found_in_list = False
			
		if not(found_in_list):
			raise Exception("ERROR: The lead SNP wasn't found in SNP ids!\n" \
				+ "Lead SNP: " + lead_snp_like_object_with_z_score.snp_id + "\n" \
				+ "SNP_ids:" + "\n" \
				+ json.dumps(SNP_ids) \
			)
	
	snps_with_approximated_z_scores = compute_approximated_zscores_for_snps_from_multiple_lead_snps(
		ld_correlation_matrix = r2_array,
		SNP_ids               = SNP_ids,
		lead_snps             = lead_snp_like_objects_with_z_scores,
	)
	
	# For every snp, check, if there is a gwas version of it with a measured 
	# pvalue.
	#
	# If there is, replace the snp with the gwas_snp.
	#
	for index, snp_with_approximated_gwas_zscore in enumerate(snps_with_approximated_z_scores):
		
		for gwas_snp in z_scores_for_lead_snp_like_objects:
			
			if gwas_snp.snp.rsID == snp_with_approximated_gwas_zscore.rsID:
				
				snps_with_approximated_z_scores[index] = gwas_snp
	
	return snps_with_approximated_z_scores, r2_array, SNP_ids

def compute_approximated_zscores_for_snps_from_multiple_lead_snps(
		ld_correlation_matrix, 
		lead_snps,
		SNP_ids
	):
	
	logger = logging.getLogger(__name__)
	
	number_of_ld_snps = len(SNP_ids)
	
	import numpy
	maximum_zscores = numpy.zeros(number_of_ld_snps)
	
	for lead_snp in lead_snps:
		approximated_zscores_for_snps_from_lead_snp = compute_approximated_zscores_for_snps_from_lead_snp(
			ld_correlation_matrix = ld_correlation_matrix,
			SNP_ids               = SNP_ids,
			lead_snp              = lead_snp,
		)
		
		for x in range(number_of_ld_snps):
			if abs(approximated_zscores_for_snps_from_lead_snp[x]) > abs(maximum_zscores[x]):
				maximum_zscores[x] = approximated_zscores_for_snps_from_lead_snp[x]
	
	from postgap.DataModel import SNP
	
	snps_with_approximated_zscores = []
	
	for index, maximum_zscore in enumerate(maximum_zscores):

		snps_with_approximated_zscores.append(
			SNP(
				rsID = SNP_ids[index],
				chrom = None,
				pos = None,
				approximated_zscore = maximum_zscore,
			)
		)
	return snps_with_approximated_zscores

def compute_approximated_zscores_for_snps_from_lead_snp(
		ld_correlation_matrix, 
		lead_snp,
		SNP_ids
	):
	
	try:
		index_of_gwas_snp = SNP_ids.index(lead_snp.snp.rsID)
	except ValueError:
		error_message = "The lead SNP wasn't found in SNP ids!\n" \
			+ "ERROR: lead SNP: %s" % lead_snp.snp_id + "\n" \
			+ "ERROR: SNP_ids:" \
			+ json.dumps(SNP_ids)
		raise IntegrityCheckException(error_message)
	
	return compute_approximated_effect_size_from_ld(
		ld_correlation_matrix = ld_correlation_matrix,
		index_of_gwas_snp     = index_of_gwas_snp, 
		lead_snp_effect_size  = lead_snp.z_score
	)

def compute_approximated_effect_size_from_ld(ld_correlation_matrix, index_of_gwas_snp, lead_snp_effect_size):
	'''
		Email Verena from 01/06/2017:
		
		```
			following

			https://www.ncbi.nlm.nih.gov/pubmed/26189819

			and the attached formulas this should be the python code

			#INPUT
			####
			#m integer number of SNPs in region
			#z_top float effect size of lead SNP
			#pos_top integer position in region of lead SNP (1,...,m)
			#LD correlation matrix (m times m)

			#OUTPUT
			#####
			#z_c array effect size vector of region of interest of length (m) interpolated from the lead SNP and regional LD

			import numpy

			z_causal = numpy.zeros(m)
			z_causal[pos_top] = z_top

			z_c=numpy.dot(LD,z_causal)
		```
	'''
	# number of SNPs in region
	m = len(ld_correlation_matrix)
	
	# Effect size of lead SNP
	z_top = lead_snp_effect_size
	
	# Position in region of lead SNP (1,..,m)
	pos_top = index_of_gwas_snp
	
	import numpy

	z_causal = numpy.zeros(m)
	z_causal[pos_top] = z_top

	z_c = numpy.dot(ld_correlation_matrix, z_causal)
	
	return z_c

def stringify_vector(vector):
	return json.dumps(vector, indent=4)
	
def stringify_matrix(matrix):
	
	import numpy as np
	# https://docs.scipy.org/doc/numpy/reference/generated/numpy.set_printoptions.html
	np.set_printoptions(precision=3, suppress=True, linewidth=150)
	return str(np.matrix(matrix))
