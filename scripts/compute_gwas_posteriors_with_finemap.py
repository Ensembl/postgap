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

python scripts/compute_gwas_posteriors_with_finemap.py \
    --gwas_cluster_file      /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir/file_based_gwas//Alzheimers/gwas_clusters/gwas_cluster_with_10_snps_around_rs11055612.pvalues_from_file.pickle \
    --output_posteriors_file /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir/file_based_gwas//Alzheimers/gwas_clusters/gwas_cluster_with_10_snps_around_rs11055612.pvalues_from_file.posteriors.pickle


python scripts/load_pvalues_for_gwas_ld_snps.py \
    --diseases \
    --gwas_cluster_file_in  /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir/EFO_0000203/gwas_clusters/gwas_cluster_with_803_snps_around_rs73071352.pickle \
    --gwas_cluster_file_out /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir/EFO_0000203/gwas_clusters/gwas_cluster_with_803_snps_around_rs73071352.pvalues_from_file.pickle

python scripts/compute_gwas_posteriors_with_finemap.py \
    --gwas_cluster_file      /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir/EFO_0000203/gwas_clusters/gwas_cluster_with_803_snps_around_rs73071352.pvalues_from_file.pickle \
    --output_posteriors_file /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir/EFO_0000203/gwas_clusters/gwas_cluster_with_803_snps_around_rs73071352.posteriors.pickle

python scripts/compute_eqtl_posteriors_with_finemap.py \
    --gwas_cluster_file      /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir/EFO_0000203/gwas_clusters/gwas_cluster_with_803_snps_around_rs73071352.pickle \
    --eqtl_cluster_file      /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir/file_based_gwas//Alzheimers/gwas_clusters/gwas_cluster_with_104_snps_around_rs10948363/cis_regulatory_evidence_from_eqtl_84_snps_linked_to_ENSG00000164393_in_Esophagus_Mucosa.pickle \
    --output_posteriors_file /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir/file_based_gwas//Alzheimers/gwas_clusters/gwas_cluster_with_104_snps_around_rs10948363/cis_regulatory_evidence_from_eqtl_84_snps_linked_to_ENSG00000164393_in_Esophagus_Mucosa.joint_posteriors.pickle




python scripts/compute_gwas_posteriors_with_finemap.py \
    --gwas_cluster_file      /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir/EFO_0000203/gwas_clusters/gwas_cluster_with_155_snps_around_rs138740.pickle \
    --output_posteriors_file /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir/EFO_0000203/gwas_clusters/gwas_cluster_with_155_snps_around_rs138740.posteriors.pickle

python scripts/compute_eqtl_posteriors_with_finemap.py \
    --gwas_cluster_file      /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir/EFO_0000203/gwas_clusters/gwas_cluster_with_155_snps_around_rs138740.pickle \
    --eqtl_cluster_file      /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir/EFO_0000203/gwas_clusters/gwas_cluster_with_155_snps_around_rs138740/cis_regulatory_evidence_from_eqtl_131_snps_linked_to_ENSG00000100284_in_Muscle_Skeletal.pickle \
    --output_posteriors_file /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir/EFO_0000203/gwas_clusters/gwas_cluster_with_155_snps_around_rs138740/cis_regulatory_evidence_from_eqtl_131_snps_linked_to_ENSG00000100284_in_Muscle_Skeletal.posteriors.pickle

python scripts/compute_joint_posteriors.py \
    --gwas_posterior          /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir/EFO_0000203/gwas_clusters/gwas_cluster_with_155_snps_around_rs138740.posteriors.pickle \
    --eqtl_posterior          /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir/EFO_0000203/gwas_clusters/gwas_cluster_with_155_snps_around_rs138740/cis_regulatory_evidence_from_eqtl_131_snps_linked_to_ENSG00000100284_in_Muscle_Skeletal.posteriors.pickle \
    --output_joint_posteriors /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir/EFO_0000203/gwas_clusters/gwas_cluster_with_155_snps_around_rs138740/cis_regulatory_evidence_from_eqtl_131_snps_linked_to_ENSG00000100284_in_Muscle_Skeletal.joint_posteriors.pickle



	'''
	
	import postgap.Globals
	postgap.Globals.DATABASES_DIR = '/nfs/nobackup/ensembl/mnuhn/postgap/databases/'

	parser = argparse.ArgumentParser(description='Run finemap')
	parser.add_argument('--gwas_cluster_file')
	parser.add_argument('--output_posteriors_file')
	
	options = parser.parse_args()
	
	logger.info( "gwas_cluster_file      = " + options.gwas_cluster_file      )
	logger.info( "output_posteriors_file = " + options.output_posteriors_file )
	
	import pickle
	
	pickle_fh    = open(options.gwas_cluster_file, 'rb')
	gwas_cluster = pickle.load(pickle_fh)
	
	from postgap.Summarisers import summarise
	logger.info("Gwas cluster:\n" + summarise(gwas_cluster))
	
	#
	# finemap/EFO_0000203/gwas_clusters/gwas_cluster_0.posteriors.pickle -> gwas_cluster_0.posteriors
	#
	import os
	# base = gwas_cluster_0.posteriors.pickle
	base = os.path.basename(options.gwas_cluster_file)
	# temp = gwas_cluster_0.posteriors
	temp = os.path.splitext(base)[0]
	# title = gwas_cluster_0
	title = os.path.splitext(temp)[0]
	
	logger.info("Cluster title: " + title)

	from postgap.FinemapIntegration.GWAS_Cluster import compute_finemap_posteriors
	gwas_posteriors = compute_finemap_posteriors(
		lead_snps    = gwas_cluster.gwas_snps,
		ld_snps      = gwas_cluster.ld_snps,
		cluster_name = title
	)
	
	if gwas_posteriors is None:
		from postgap.FinemapIntegration.GWAS_Cluster import FinemapFailedException
		raise FinemapFailedException("Got no finemap results!")

	logging.info( "Gwas posteriors have been computed:\n" + summarise(gwas_posteriors) )
	
	import os
	output_posteriors_file = options.output_posteriors_file
	output_directory = os.path.dirname(output_posteriors_file)
	
	if not os.path.exists(output_directory):
		os.makedirs(output_directory)

	pickle_fh = open(output_posteriors_file, 'w')
	pickle.dump(gwas_posteriors, pickle_fh)
	pickle_fh.close
	
	logging.info("Posteriors have been written to " + output_posteriors_file)
	logging.info("All done.");

if __name__ == "__main__":
	main()

