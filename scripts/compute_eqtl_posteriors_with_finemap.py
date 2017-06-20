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
import sys
import argparse
import collections
import json
import sqlite3
import logging
import logging.config

logging.config.fileConfig('configuration/logging.conf')
logger = logging.getLogger(__name__)

def main():
	'''

python scripts/compute_eqtl_posteriors_with_finemap.py \
    --eqtl_cluster_file      finemap/EFO_0000203/eqtl_clusters/eqtl_snps_linked_to_ENSG00000168038_in_Whole_Blood.pickle \
    --output_posteriors_file finemap/EFO_0000203/eqtl_clusters/eqtl_snps_linked_to_ENSG00000168038_in_Whole_Blood.posteriors.pickle

	'''
	
	import postgap.Globals
	postgap.Globals.DATABASES_DIR = '/nfs/nobackup/ensembl/mnuhn/postgap/databases/'

	parser = argparse.ArgumentParser(description='Run finemap')
	parser.add_argument('--eqtl_cluster_file')
	parser.add_argument('--output_posteriors_file')
	
	options = parser.parse_args()

	logger.info( "eqtl_cluster_file      = " + options.eqtl_cluster_file      )
	logger.info( "output_posteriors_file = " + options.output_posteriors_file )
	
	import pickle
	
	pickle_fh       = open(options.eqtl_cluster_file, 'rb')
	eqtl_cluster    = pickle.load(pickle_fh)
	
	eqtl_posteriors = process_eqtl_cluster(eqtl_cluster)
	
	import os
	output_posteriors_file = options.output_posteriors_file
	output_directory = os.path.dirname(output_posteriors_file)
	
	if not os.path.exists(output_directory):
		os.makedirs(output_directory)

	pickle_fh = open(output_posteriors_file, 'w')
	pickle.dump(eqtl_posteriors, pickle_fh)
	pickle_fh.close
	
	logging.info("Posteriors have been written to " + output_posteriors_file)

	logging.info("All done.");

def process_eqtl_cluster(eqtl_cluster):
	
	import pprint
	pp = pprint.PrettyPrinter(indent=4, width=80)
	pp.pprint(eqtl_cluster)
	
	snps_correlated_to_current_tissue_and_gene = []

	for cisregulatory_evidence in eqtl_cluster:
		snps_correlated_to_current_tissue_and_gene.append(cisregulatory_evidence.snp)
	
	# Finemap needs zscores.
	#
	# In order to compute zscores the betas are needed in addition 
	# to the pvalues.
	#
	# The betas are fetched from the rest server.
	#
	zscore_vector = []
	
	snp_hash = dict( (snp.rsID, snp) for snp in snps_correlated_to_current_tissue_and_gene)
	
	from postgap.Cisreg import GTEx
	GTEx = GTEx()
	#
	# At this point all cisregulatory evidence should be GTEx and from the 
	# same snp and tissue.
	#
	gene   = eqtl_cluster[0].gene
	tissue = eqtl_cluster[0].tissue
	cisregulatory_evidence_with_betas = GTEx.gene_tissue_betas(gene, tissue, snp_hash)
	
	snp_to_beta = dict(
		(cisregulatory_evidence_with_beta.snp.rsID, cisregulatory_evidence_with_beta.beta) 
			for cisregulatory_evidence_with_beta in cisregulatory_evidence_with_betas
	)
	
	# The betas are in snp_to_beta now.
	
	print "snp_to_beta:"
	pp.pprint(snp_to_beta)
	
	# Compute vector of zscores.
	#
	
	from methods.GWAS_SNP import compute_z_score_from_pvalue_and_sign, sign
	
	for cisregulatory_evidence in eqtl_cluster:
		pvalue = cisregulatory_evidence.pvalue
		zscore = compute_z_score_from_pvalue_and_sign(
			pvalue, 
			sign(snp_to_beta[cisregulatory_evidence.snp.rsID])
		)
		zscore_vector.append(zscore)
	
	print "zscore_vector:"
	pp.pprint(zscore_vector)
	
	print "snps_correlated_to_current_tissue_and_gene:"
	pp.pprint(snps_correlated_to_current_tissue_and_gene)
	
	from postgap.LD import get_pairwise_ld
	(SNP_ids, r2_array) = get_pairwise_ld(snps_correlated_to_current_tissue_and_gene)
	
	logger.info("SNPs:")
	logger.info("\n" + stringify_vector(SNP_ids))
	
	logger.info("zscores:")
	logger.info("\n" + stringify_vector(zscore_vector))
	
	logger.info("Matrix of pairwise linkage disequilibria:")
	logger.info("\n" + stringify_matrix(r2_array))

	logger.info("Running finemap")
	
	kstart = 1
	kmax   = 1
	max_iter = "Not used when kstart == kmax"

	import finemap.stochastic_search as sss
	finemap_posteriors = sss.finemap(
		z_scores   = zscore_vector,
		cov_matrix = r2_array,
		n          = len(r2_array),
		kstart     = kstart,
		kmax       = kmax,
		max_iter   = max_iter,
		prior      = "independence"
	)
	logger.info("Done running finemap")
	
	if finemap_posteriors is None:
		raise Exception ("ERROR: Didn't get posteriors!")
	
	return finemap_posteriors

def stringify_vector(vector):
	import json
	return json.dumps(vector, indent=4)
	
def stringify_matrix(matrix):
	
	import numpy as np
	# https://docs.scipy.org/doc/numpy/reference/generated/numpy.set_printoptions.html
	np.set_printoptions(precision=3, suppress=True, linewidth=150)
	return str(np.matrix(matrix))


if __name__ == "__main__":
	main()

