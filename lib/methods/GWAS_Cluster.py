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
snp_with_zscore = collections.namedtuple(
	'snp_with_zscore', 
	[
		'snp_id',
		'z_score'
	]
)

class ZScoreComputationException(Exception):
	pass

class FinemapFailedException(Exception):
	pass

class IntegrityCheckException(Exception):
	pass


def GWAS_Clusters_ok(gwas_clusters):
	for gwas_cluster in gwas_clusters:
		if not(ld_snps_contain_gwas_snps(gwas_cluster)):
			return False, "Gwas snps not among ld snps"
	return True, ""

def compute_gwas_cluster_with_finemap_posteriors(gwas_cluster, cluster_name="Unnamed cluster"):
	
	try:
		finemap_posteriors = compute_finemap_posteriors(gwas_cluster, cluster_name = cluster_name)
		
		if finemap_posteriors is None:
			raise FinemapFailedException("Got no finemap results!")

	except ZScoreComputationException, z:
		logger = logging.getLogger(__name__)
		logger.warning(str(z))
		logger.warning("Skipping this cluster.")
		return
	
	return GWAS_Cluster(
		gwas_snps          = gwas_cluster.gwas_snps,
		ld_snps            = gwas_cluster.ld_snps,
		finemap_posteriors = finemap_posteriors
	)

def compute_finemap_posteriors(gwas_cluster, cluster_name="Unnamed cluster"):
	
	ld_snps   = gwas_cluster.ld_snps
	gwas_snps = gwas_cluster.gwas_snps
	
	(approximated_gwas_zscores, r2_array, SNP_ids) = compute_approximated_gwas_zscores(
		gwas_snps = gwas_snps,
		ld_snps   = ld_snps
	)
	logger = logging.getLogger(__name__)
	logger.info("Running finemap")
	
	kstart = 1
	kmax   = 1
	max_iter = "Not used when kstart == kmax"

	import finemap.stochastic_search as sss
	finemap_posteriors = sss.finemap(
		labels       = SNP_ids,
		z_scores     = approximated_gwas_zscores,
		cov_matrix   = r2_array,
		n            = len(r2_array),
		kstart       = kstart,
		kmax         = kmax,
		max_iter     = max_iter,
		prior        = "independence",
		sample_label = cluster_name
	)
	logger.info("Done running finemap")
	
	if finemap_posteriors is None:
		raise Exception ("ERROR: Didn't get posteriors!")

	return finemap_posteriors

def compute_gwas_snps_with_z_scores(gwas_snps):
	
	from methods.GWAS_SNP import compute_z_score_for_gwas_snp
	
	gwas_snps_with_z_scores = filter (
		lambda X: X.z_score is not None, [ 
			snp_with_zscore(
				snp_id  = gwas_snp.snp.rsID,
				z_score = compute_z_score_for_gwas_snp(gwas_snp)
			)
				for gwas_snp in gwas_snps 
		]
	)
	gwas_cluster_size = len(gwas_snps)
	
	logger = logging.getLogger(__name__)
	
	if len(gwas_snps_with_z_scores) == 0:
		raise ZScoreComputationException("Couldn't compute a z score for any of the gwas snps.")
	
	return gwas_snps_with_z_scores

def compute_approximated_gwas_zscores(gwas_snps, ld_snps):
	# LD measures the deviation from the expectation of non-association.
	(SNP_ids, r2_array) = postgap.LD.get_pairwise_ld(ld_snps)

	logger = logging.getLogger(__name__)
	
	ld_cluster_size = len(ld_snps)

	if ld_cluster_size==0:
		raise Exception("The cluster size was zero!")

	logger.info("SNPs:")
	logger.info("\n" + stringify_vector(SNP_ids))

	logger.info("Matrix of pairwise linkage disequilibria:")
	logger.info("\n" + stringify_matrix(r2_array))
	
	gwas_snps_with_z_scores = compute_gwas_snps_with_z_scores(gwas_snps)
	
	# This loop is just a healthcheck
	for gwas_snps_with_z_score in gwas_snps_with_z_scores:
		
		found_in_list = None
		try:
			SNP_ids.index(gwas_snps_with_z_score.snp_id)
			found_in_list = True
		except ValueError:
			found_in_list = False
			
		if not(found_in_list):
			raise Error("ERROR: The lead SNP wasn't found in SNP ids!\n" \
				+ "lead SNP: " + gwas_snps_with_z_score.snp_id + "\n" \
				+ "SNP_ids:" + "\n" \
				+ json.dumps(SNP_ids) \
			)
	
	from methods.GWAS_Cluster import compute_approximated_zscores_for_snps_from_multiple_lead_snps
	
	approximated_gwas_zscores = compute_approximated_zscores_for_snps_from_multiple_lead_snps(
		ld_correlation_matrix = r2_array,
		SNP_ids               = SNP_ids,
		lead_snps             = gwas_snps_with_z_scores,
	)
	
	return approximated_gwas_zscores, r2_array, SNP_ids

def ld_snps_contain_gwas_snps(gwas_cluster):
	
	gwas_snps = gwas_cluster.gwas_snps
	ld_snps   = gwas_cluster.ld_snps
	
	for gwas_snp in gwas_snps:
		
		snp  = gwas_snp.snp
		rsID = snp.rsID
		
		if not( any([ ld_snp.rsID == rsID for ld_snp in ld_snps ]) ):
			return False
	return True

def compute_approximated_zscores_for_snps_from_multiple_lead_snps(
		ld_correlation_matrix, 
		lead_snps,
		SNP_ids
	):

	number_of_ld_snps = len(SNP_ids)
	
	import numpy
	maximum_zscores = numpy.zeros(number_of_ld_snps)
	
	for lead_snp in lead_snps:
		approximated_zscores_for_snps_from_lead_snp = compute_approximated_zscores_for_snps_from_lead_snp(
			ld_correlation_matrix = ld_correlation_matrix,
			SNP_ids               = SNP_ids,
			lead_snp              = lead_snp,
		)

		#print "approximated_zscores_for_snps_from_lead_snp = "
		#pprint.pprint(approximated_zscores_for_snps_from_lead_snp)
		
		for x in range(number_of_ld_snps):
			if abs(approximated_zscores_for_snps_from_lead_snp[x]) > abs(maximum_zscores[x]):
				maximum_zscores[x] = approximated_zscores_for_snps_from_lead_snp[x]

	#print "final zscores = "
	#pprint.pprint(maximum_zscores)
	
	return maximum_zscores

def compute_approximated_zscores_for_snps_from_lead_snp(
		ld_correlation_matrix, 
		lead_snp,
		SNP_ids
	):
	
	try:
		index_of_gwas_snp = SNP_ids.index(lead_snp.snp_id)
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

#print compute_z_score_from_pvalue_and_odds_ratio_or_beta_coefficient(0.0001, -1)
#print compute_z_score_from_pvalue_and_odds_ratio_or_beta_coefficient(0.0001, beta_coefficient=-1)
#print compute_z_score_from_pvalue_and_odds_ratio_or_beta_coefficient(0.001, beta_coefficient=-1)
#print compute_z_score_from_pvalue_and_odds_ratio_or_beta_coefficient(0.001)

#sys.exit(0)

def stringify_vector(vector):
	return json.dumps(vector, indent=4)
	
def stringify_matrix(matrix):
	
	import numpy as np
	# https://docs.scipy.org/doc/numpy/reference/generated/numpy.set_printoptions.html
	np.set_printoptions(precision=3, suppress=True, linewidth=150)
	return str(np.matrix(matrix))
