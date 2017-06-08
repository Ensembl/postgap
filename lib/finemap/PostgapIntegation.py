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
import pprint, pickle

class ZScoreComputationException(Exception):
	pass

class FinemapFailedException(Exception):
	pass

class IntegrityCheckException(Exception):
	pass

def process_GWAS_Clusters(gwas_clusters):
	
	assert(GWAS_Clusters_ok(gwas_clusters))
	
	gwas_clusters_with_posteriors = [ process_GWAS_Cluster(gwas_cluster) for gwas_cluster in gwas_clusters ]
	return gwas_clusters_with_posteriors

def GWAS_Clusters_ok(gwas_clusters):
	for gwas_cluster in gwas_clusters:
		if not(ld_snps_contain_gwas_snps(gwas_cluster)):
			return False, "Gwas snps not among ld snps"
	return True, ""

def ld_snps_contain_gwas_snps(gwas_cluster):
	
	gwas_snps = gwas_cluster.gwas_snps
	ld_snps   = gwas_cluster.ld_snps
	
	for gwas_snp in gwas_snps:
		
		snp  = gwas_snp.snp
		rsID = snp.rsID
		
		if not( any([ ld_snp.rsID == rsID for ld_snp in ld_snps ]) ):
			return False
	return True


def process_GWAS_Cluster(gwas_cluster):
	
	try:
		finemap_posteriors = process_ld_snps(
			ld_snps   = gwas_cluster.ld_snps,
			gwas_snps = gwas_cluster.gwas_snps
		)
		if finemap_posteriors is None:
			raise FinemapFailedException("Got no finemap results!")

	except ZScoreComputationException, z:
		print "WARNING: " + str(z)
		print "WARNING: Skipping this cluster."
		return
	
	return GWAS_Cluster(
		gwas_snps          = gwas_cluster.gwas_snps,
		ld_snps            = gwas_cluster.ld_snps,
		finemap_posteriors = finemap_posteriors
	)

def compute_z_score_from_pvalue_and_sign(pvalue, sign):

	from scipy.stats import norm
	
	if sign is None:
		return None

	# norm.ppf means:
	# Percent point function (inverse of cdf) at q of the given RV.
	#
	# See:
	# https://docs.scipy.org/doc/scipy-0.14.0/reference/stats.html
	
	z_score = - norm.ppf(pvalue/2) * sign
	
	return z_score

def compute_z_score_sign(odds_ratio, beta_coefficient):
	'''
		where the gwas_sign comes from

		beta > 0 gwas_sign = +1
		beta < 0 gwas_sign = -1

		or > 1 gwas_sign = +1
		or < 1 gwas_sign = -1
	'''
	
	if beta_coefficient>0:
		return 1
	
	# None<0 is true, so must test for that separately.
	if beta_coefficient is not None and beta_coefficient<0:
		return -1
	
	if odds_ratio>0:
		return 1
	
	# None<0 is true, so must test for that separately.
	if odds_ratio is not None and odds_ratio<0:
		return -1
	
	return None;

def compute_z_score_from_pvalue_and_odds_ratio_or_beta_coefficient(pvalue, odds_ratio=None, beta_coefficient=None):
	
	return compute_z_score_from_pvalue_and_sign(
		pvalue,
		compute_z_score_sign(odds_ratio, beta_coefficient)
	)

def compute_z_score_for_gwas_snp(gwas_snp):
	
	return compute_z_score_from_pvalue_and_odds_ratio_or_beta_coefficient(
		pvalue           = gwas_snp.pvalue,
		odds_ratio       = gwas_snp.odds_ratio,
		beta_coefficient = gwas_snp.beta_coefficient,
	)

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

def print_vector(vector):
	print json.dumps(vector, indent=4)
	
def print_matrix(matrix):
	
	import numpy as np
	# https://docs.scipy.org/doc/numpy/reference/generated/numpy.set_printoptions.html
	np.set_printoptions(precision=3, suppress=True, linewidth=150)
	print(np.matrix(matrix))

def process_ld_snps(ld_snps, gwas_snps):

	if len(gwas_snps)!=1:
		print "WARNING: Got a cluster with more than one (%i) gwas SNP, but only looking at the first!" % len(gwas_snps)

	ld_cluster_size = len(ld_snps)

	if ld_cluster_size==1:
		print "The cluster size was one!"
	if ld_cluster_size==0:
		print "ERROR: The cluster size was zero!"

	print "The cluster has %i ld members.\n" % ld_cluster_size

	#gwas_snp = gwas_snps[0]
	#print "gwas snp: %s" % gwas_snp.snp.rsID
	#print "  - pvalue:                     %s" % gwas_snp.pvalue
	#print "  - odds_ratio:                 %s" % gwas_snp.odds_ratio
	#print "  - beta_coefficient:           %s" % gwas_snp.beta_coefficient
	#print "  - beta_coefficient_unit:      %s" % gwas_snp.beta_coefficient_unit
	#print "  - beta_coefficient_direction: %s" % gwas_snp.beta_coefficient_direction
		
	import collections
	snp_with_zscore = collections.namedtuple(
		'snp_with_zscore', 
		[
			'snp_id',
			'z_score'
		]
	)
	
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
	
	print "The cluster has %i gwas members.\n" % gwas_cluster_size
	print "The cluster has %i gwas members with z scores.\n" % len(gwas_snps_with_z_scores)

	if len(gwas_snps_with_z_scores) == 0:
		raise ZScoreComputationException("Couldn't compute a z score for any of the gwas snps.")

	print "Location of the SNPs:"
	
	snps_from_gwas_snps = [ gwas_snp.snp for gwas_snp in gwas_snps ]
	
	print "---- Gwas SNPs ----"
	for snp in snps_from_gwas_snps:
		print "    %s    %s    %s" % (snp.rsID, snp.chrom, snp.pos)

	# The GWAS SNP is always one of the LD SNPs.
	print "---- LD SNPs ----"
	for snp in ld_snps:
		print "    %s    %s    %s" % (snp.rsID, snp.chrom, snp.pos)

	# LD measures the deviation from the expectation of non-association.
	(SNP_ids, r2_array) = postgap.LD.get_pairwise_ld(ld_snps)

	print "SNPs:\n"
	print_vector(SNP_ids)

	print "Matrix of pairwise linkage disequilibria:\n"
	#set_diagonal(r2_array)
	print_matrix(r2_array)
	
	for gwas_snps_with_z_score in gwas_snps_with_z_scores:
		
		found_in_list = None
		try:
			SNP_ids.index(gwas_snps_with_z_score.snp_id)
			found_in_list = True
		except ValueError:
			found_in_list = False
			
		if not(found_in_list):
			print "ERROR: The lead SNP wasn't found in SNP ids!"
			print "ERROR: lead SNP: %s" % gwas_snps_with_z_score.snp_id
			print "ERROR: SNP_ids:"
			print json.dumps(SNP_ids)
			sys.exit(1)
	
	approximated_gwas_zscore = compute_approximated_zscores_for_snps_from_multiple_lead_snps(
		ld_correlation_matrix = r2_array,
		SNP_ids               = SNP_ids,
		lead_snps             = gwas_snps_with_z_scores,
	)
	print "Running finemap"
	
	kstart = 3
	kmax   = 3
	max_iter = "Not used when kstart == kmax"

	import finemap.stochastic_search as sss
	finemap_posteriors = sss.finemap(
		z_scores   = approximated_gwas_zscore,
		cov_matrix = r2_array,
		n          = len(r2_array),
		kstart     = kstart,
		kmax       = kmax,
		max_iter   = max_iter,
		prior      = "independence"
	)
	#pprint.pprint(finemap_posteriors)
	print "Done running finemap"
	
	if finemap_posteriors is None:
		raise Exception ("ERROR: Didn't get posteriors!")

	return finemap_posteriors

