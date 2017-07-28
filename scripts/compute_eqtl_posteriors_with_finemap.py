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
    --gwas_cluster_file      /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir_for_diabetes5/diabetes/gwas_clusters/gwas_cluster_around_snp_rs10510110.pickle \
    --eqtl_cluster_file      /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir_for_diabetes5/diabetes/gwas_clusters/gwas_cluster_around_snp_rs10510110/cis_regulatory_evidence_from_eqtl_snp_rs11200594_linked_to_ENSG00000107679_in_Whole_Blood.pickle \
    --output_posteriors_file /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir_for_diabetes5/diabetes/gwas_clusters/gwas_cluster_around_snp_rs10510110/cis_regulatory_evidence_from_eqtl_snp_rs11200594_linked_to_ENSG00000107679_in_Whole_Blood.posteriors.pickle

python scripts/stringify_pickle_object.py --pickle_file /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir_for_diabetes5/diabetes/gwas_clusters/gwas_cluster_around_snp_rs10510110/cis_regulatory_evidence_from_eqtl_snp_rs2034245_linked_to_ENSG00000155542_in_Whole_Blood.pickle

	'''
	
	import postgap.Globals
	postgap.Globals.DATABASES_DIR = '/nfs/nobackup/ensembl/mnuhn/postgap/databases/'

	parser = argparse.ArgumentParser(description='Run finemap')
	parser.add_argument('--gwas_cluster_file')
	parser.add_argument('--eqtl_cluster_file')
	parser.add_argument('--output_posteriors_file')
	
	options = parser.parse_args()

	logging.info( "gwas_cluster_file      = " + options.gwas_cluster_file      )
	logging.info( "eqtl_cluster_file      = " + options.eqtl_cluster_file      )
	logging.info( "output_posteriors_file = " + options.output_posteriors_file )
	
	import pickle
	
	pickle_fh       = open(options.gwas_cluster_file, 'rb')
	gwas_cluster    = pickle.load(pickle_fh)
	
	assert type(gwas_cluster) is postgap.DataModel.GWAS_Cluster, "gwas_cluster is a GWAS_Cluster"

	#pickle_fh       = open(options.eqtl_cluster_file, 'rb')
	#eqtl_cluster    = pickle.load(pickle_fh)
	
	eqtl_clusters = []
	
	with open(options.eqtl_cluster_file, 'rb') as pickle_fh:
		try:
			while True:
				eqtl_cluster = pickle.load(pickle_fh)
				assert type(eqtl_cluster) is postgap.DataModel.Cisregulatory_Evidence, "eqtl_cluster is a Cisregulatory_Evidence"
				eqtl_clusters.append(eqtl_cluster)
		except EOFError:
			logger.info("Got " + str(len(eqtl_clusters)) + " eqtl clusters.")
	
	for eqtl_cluster in eqtl_clusters:
		from pprint import pformat
		logger.info(pformat(eqtl_cluster))
	
	logger.info(pformat(gwas_cluster))
	
	assert_all_eqtl_snps_are_in_gwas_cluster(
		eqtl_clusters = eqtl_clusters,
		gwas_cluster = gwas_cluster
	)

	#
	# finemap/EFO_0000203/gwas_clusters/gwas_cluster_0.posteriors.pickle -> gwas_cluster_0.posteriors
	#
	import os
	# base = gwas_cluster_0.posteriors.pickle
	base = os.path.basename(options.eqtl_cluster_file)
	# temp = gwas_cluster_0.posteriors
	temp = os.path.splitext(base)[0]
	# title = gwas_cluster_0
	title = os.path.splitext(temp)[0]

	eqtl_posteriors = process_eqtl_cluster(
		cisregulatory_evidence = eqtl_cluster, 
		gwas_cluster = gwas_cluster,
		cluster_name = title
	)
	
	logging.info( "eQTL posteriors have been computed:\n" + str(eqtl_posteriors) )
	
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

def assert_all_eqtl_snps_are_in_gwas_cluster(gwas_cluster, eqtl_clusters):
	
	all_snp_ids_in_gwas_cluster = dict()
	
	for gwas_snp in gwas_cluster.gwas_snps:
		all_snp_ids_in_gwas_cluster[gwas_snp.snp.rsID] = 1
	
	for ld_snp in gwas_cluster.ld_snps:
		all_snp_ids_in_gwas_cluster[ld_snp.rsID] = 1
	
	for eqtl_cluster in eqtl_clusters:

		if not(eqtl_cluster.snp.rsID in all_snp_ids_in_gwas_cluster):
			raise Exception("Found id in eqtl cluster that is not in the gwas cluster!")
		
		#print "Ok " + eqtl_cluster.snp.rsID	

def process_eqtl_cluster(cisregulatory_evidence, gwas_cluster, cluster_name = "Unnamed cluster"):
	
	import pprint
	pp = pprint.PrettyPrinter(indent=4, width=80)
	pp.pprint(cisregulatory_evidence)
	
	ld_snps   = gwas_cluster.ld_snps
	gwas_snps = gwas_cluster.gwas_snps
	
	snps_from_gwas_snps = [ gwas_snp.snp for gwas_snp in gwas_snps ]
	
	from methods.GWAS_Cluster import compute_approximated_gwas_zscores
	(approximated_gwas_zscores, r2_array, SNP_ids) = compute_approximated_gwas_zscores(
		gwas_snps = [ cisregulatory_evidence ],
		ld_snps   = ld_snps + snps_from_gwas_snps
	)
	
	logging.info("SNPs:")
	logging.info("\n" + stringify_vector(SNP_ids))
	
	logging.info("zscores " + str(type(approximated_gwas_zscores)) + ":")
	logging.info("\n" + stringify_vector(approximated_gwas_zscores))
	
	logging.info("Matrix of pairwise linkage disequilibria:")
	logging.info("\n" + stringify_matrix(r2_array))

	logging.info("Cluster name: " + cluster_name)

	logging.info("Running finemap")

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
	logging.info("Done running finemap")
	
	if finemap_posteriors is None:
		raise Exception ("ERROR: Didn't get posteriors!")
	
	return finemap_posteriors

def stringify_vector(vector):
	#import json
	#return json.dumps(vector, indent=4)
	from pprint import pformat
	return pformat(vector)
	
def stringify_matrix(matrix):
	
	import numpy as np
	# https://docs.scipy.org/doc/numpy/reference/generated/numpy.set_printoptions.html
	np.set_printoptions(precision=3, suppress=True, linewidth=150)
	return str(np.matrix(matrix))


if __name__ == "__main__":
	main()

