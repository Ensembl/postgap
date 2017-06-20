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
import logging
import logging.config

logging.config.fileConfig('configuration/logging.conf')
logger = logging.getLogger(__name__)

def main():
	'''
	
python scripts/compute_joint_posteriors.py \
    --gwas_posterior          finemap/EFO_0000203/gwas_clusters_posteriors/gwas_cluster_0.posteriors.pickle \
    --eqtl_posterior          finemap/EFO_0000203/eqtl_clusters_posteriors/eqtl_snps_linked_to_ENSG00000168038_in_Whole_Blood.posteriors.pickle \
    --output_joint_posteriors finemap/EFO_0000203/join_posteriors/eqtl_snps_linked_to_ENSG00000168038_in_Whole_Blood__gwas_cluster_0.joint_posteriors.pickle

	'''
	
	import postgap.Globals
	postgap.Globals.DATABASES_DIR = '/nfs/nobackup/ensembl/mnuhn/postgap/databases/'

	parser = argparse.ArgumentParser(description='Run finemap')
	parser.add_argument('--gwas_posterior')
	parser.add_argument('--eqtl_posterior')
	parser.add_argument('--output_joint_posteriors')
	
	options = parser.parse_args()

	gwas_posterior          = options.gwas_posterior
	eqtl_posterior          = options.eqtl_posterior
	output_joint_posteriors = options.output_joint_posteriors

	logger.info( "gwas_posterior              = " + gwas_posterior )
	logger.info( "eqtl_posterior              = " + eqtl_posterior )
	logger.info( "output_joint_posteriors     = " + output_joint_posteriors )
	
	import pickle
	
	pickle_fh       = open(eqtl_posterior, 'rb')
	eqtl_posteriors = pickle.load(pickle_fh)

	pickle_fh       = open(gwas_posterior, 'rb')
	gwas_posteriors = pickle.load(pickle_fh)
	
	logging.info("Computing joint posteriors")
	fm2d_sss = gwas_posteriors.joint_posterior(eqtl_posteriors)
	
	import os
	output_joint_posteriors_directory = os.path.dirname(output_joint_posteriors)
	
	if not os.path.exists(output_joint_posteriors_directory):
		os.makedirs(output_joint_posteriors_directory)

	pickle_fh = open(output_joint_posteriors, 'w')
	pickle.dump(fm2d_sss, pickle_fh)
	pickle_fh.close
	
	coloc_evidence = fm2d_sss[0]
	res_final      = fm2d_sss[1]
	
	configurations = res_final.configurations
	posterior      = res_final.posterior

	from pprint import pformat

	logging.info("coloc_evidence: %s"     % coloc_evidence)
	logging.info("Got %s configurations." % len(configurations))
	logging.info("configurations: "       + pformat(configurations))
	logging.info("Got %s posteriors: "    % len(posterior))
	logging.info("posterior: "            + pformat(posterior))

	logging.info("Posteriors have been written to " + output_joint_posteriors)
	logging.info("All done.");

if __name__ == "__main__":
	main()

