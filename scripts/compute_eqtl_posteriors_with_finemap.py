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
import argparse
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


python scripts/compute_eqtl_posteriors_with_finemap.py \
    --gwas_cluster_file      /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir/file_based_gwas//Alzheimers/gwas_clusters/gwas_cluster_with_104_snps_around_rs10948363.pickle \
    --eqtl_cluster_file      /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir/file_based_gwas//Alzheimers/gwas_clusters/gwas_cluster_with_104_snps_around_rs10948363/cis_regulatory_evidence_from_eqtl_84_snps_linked_to_ENSG00000164393_in_Esophagus_Mucosa.pickle \
    --output_posteriors_file /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir/file_based_gwas//Alzheimers/gwas_clusters/gwas_cluster_with_104_snps_around_rs10948363/cis_regulatory_evidence_from_eqtl_84_snps_linked_to_ENSG00000164393_in_Esophagus_Mucosa.joint_posteriors.pickle

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

	eqtl_clusters = []
	
	with open(options.eqtl_cluster_file, 'rb') as pickle_fh:
		try:
			while True:
				eqtl_cluster = pickle.load(pickle_fh)
				assert type(eqtl_cluster) is postgap.DataModel.Cisregulatory_Evidence, "eqtl_cluster is a Cisregulatory_Evidence"
				eqtl_clusters.append(eqtl_cluster)
		except EOFError:
			logger.info("Got " + str(len(eqtl_clusters)) + " eqtl clusters.")
	
	from postgap.Summarisers import summarise
	logger.info("The gwas cluster is:\n" + summarise(gwas_cluster))
	
	logger.info("The cisregulatory evidence is:\n" + summarise(eqtl_clusters))
	
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

	from postgap.FinemapIntegration.GWAS_Cluster import compute_finemap_posteriors
	
	ld_snps_with_cis_regulatory_evidence_inserted = []
	
	from postgap.DataModel import Cisregulatory_Evidence, GWAS_SNP, SNP
	
	for ld_snp in gwas_cluster.ld_snps:
		
		rsID = None
		
		if type(ld_snp) is SNP:
			rsID = ld_snp.rsID
		
		if type(ld_snp) is GWAS_SNP:
			rsID = ld_snp.snp.rsID
		
		if rsID is None:
			raise Exception(str(type(ld_snp)))
		
		cisregulatory_evidence = get_cisregulatory_evidence(rsID_to_find=rsID, cisregulatory_evidence_list=eqtl_clusters)
		
		if cisregulatory_evidence is None:
			ld_snps_with_cis_regulatory_evidence_inserted.append(ld_snp)
			continue
		
		if type(cisregulatory_evidence) is Cisregulatory_Evidence:
			ld_snps_with_cis_regulatory_evidence_inserted.append(cisregulatory_evidence)
			continue
			
		raise Exception(str(type(cisregulatory_evidence)))
	
	logging.info( "The ld snps with cis regulatory evidence inserted where possible:\n" + summarise(ld_snps_with_cis_regulatory_evidence_inserted) )
	
	finemap_posteriors = compute_finemap_posteriors(
		lead_snps    = eqtl_clusters,
		ld_snps      = ld_snps_with_cis_regulatory_evidence_inserted,
		cluster_name = title
	)
	
	if finemap_posteriors is None:
		raise Exception ("ERROR: Didn't get posteriors!")

	logging.info( "eQTL posteriors have been computed:\n" + str(finemap_posteriors) )
	
	import os
	output_posteriors_file = options.output_posteriors_file
	output_directory = os.path.dirname(output_posteriors_file)
	
	if not os.path.exists(output_directory):
		os.makedirs(output_directory)

	pickle_fh = open(output_posteriors_file, 'w')
	pickle.dump(finemap_posteriors, pickle_fh)
	pickle_fh.close
	
	logging.info("Posteriors have been written to " + output_posteriors_file)

	logging.info("All done.");

def get_cisregulatory_evidence(rsID_to_find, cisregulatory_evidence_list):
	
	for cisregulatory_evidence in cisregulatory_evidence_list:
		
		rsID = cisregulatory_evidence.snp.rsID
		if rsID_to_find == rsID:
			return cisregulatory_evidence
	
	return None

def assert_all_eqtl_snps_are_in_gwas_cluster(gwas_cluster, eqtl_clusters):
	
	all_snp_ids_in_gwas_cluster = dict()
	
	for gwas_snp in gwas_cluster.gwas_snps:
		all_snp_ids_in_gwas_cluster[gwas_snp.snp.rsID] = 1
	
	for ld_snp in gwas_cluster.ld_snps:
		all_snp_ids_in_gwas_cluster[ld_snp.rsID] = 1
	
	for eqtl_cluster in eqtl_clusters:

		if not(eqtl_cluster.snp.rsID in all_snp_ids_in_gwas_cluster):
			raise Exception("Found id in eqtl cluster that is not in the gwas cluster!")

if __name__ == "__main__":
	main()

