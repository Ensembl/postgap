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

import logging
import logging.config
logging.config.fileConfig('configuration/logging.conf')
logger = logging.getLogger(__name__)

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

def sign(value):
	if value>0:
		return 1
	
	# None<0 is true, so must test for that separately.
	if value is not None and value<0:
		return -1
	
	if value==0:
		return 0

	raise Exception

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

class snp_in_multiple_gwas_associations_exception(Exception):
	pass

def compute_z_score_for_snp_like_type(snp_like):
	
	from postgap.DataModel import GWAS_SNP
	if type(snp_like) is GWAS_SNP:
		return compute_z_score_for_gwas_snp(snp_like)
	
	from postgap.DataModel import Cisregulatory_Evidence
	if type(snp_like) is Cisregulatory_Evidence:
		zscore = compute_z_score_for_cisregulatory_evidence(snp_like)
		return compute_z_score_for_cisregulatory_evidence(snp_like)
	
	# Should never happen
	raise Exception

def compute_z_score_for_cisregulatory_evidence(cisregulatory_evidence):
	
	# For eqtls, this should never be true.
	risk_alleles_present_in_reference = False
	
	z_score = compute_z_score_from_pvalue_and_odds_ratio_or_beta_coefficient(
		pvalue           = cisregulatory_evidence.pvalue,
		odds_ratio       = None,
		beta_coefficient = cisregulatory_evidence.beta,
	)
	return z_score

def compute_z_score_for_gwas_snp(gwas_snp):
	
	evidence_list = gwas_snp.evidence
	
	if len(evidence_list) != 1:
		from postgap.Summarisers import summarise
		raise snp_in_multiple_gwas_associations_exception("This GWAS SNP was found via more than one GWAS Association:\n" + summarise(gwas_snp))
	
	evidence = evidence_list[0]
	
	from postgap.DataModel import GWAS_Association
	assert type(evidence) is GWAS_Association, "Object is GWAS_Association"
	
	# If the gwas risk allele is present in the reference, then the sign
	# for the zscore has to be inverted.
	risk_alleles_present_in_reference = evidence.risk_alleles_present_in_reference
	
	z_score_from_pvalue_and_odds_ratio_or_beta_coefficient = compute_z_score_from_pvalue_and_odds_ratio_or_beta_coefficient(
		pvalue           = evidence.pvalue,
		odds_ratio       = evidence.odds_ratio,
		beta_coefficient = evidence.beta_coefficient,
	)
	
	if z_score_from_pvalue_and_odds_ratio_or_beta_coefficient is None:
		return None
	
	if risk_alleles_present_in_reference:
		logger.info("The zscore will be inverted, because the risk allele is present in the reference.")
		zscore = -1 * z_score_from_pvalue_and_odds_ratio_or_beta_coefficient
		return zscore
	
	if not(risk_alleles_present_in_reference):
		logger.info("The risk allele is not present in the reference, so the zscore is not affected.")
		zscore = z_score_from_pvalue_and_odds_ratio_or_beta_coefficient
		return zscore
	
	# Should never happen
	raise Exception

