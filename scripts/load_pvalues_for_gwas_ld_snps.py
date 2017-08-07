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
from pprint import pformat

logging.config.fileConfig('configuration/logging.conf')
logger = logging.getLogger(__name__)

def main():
	'''

python scripts/load_pvalues_for_gwas_ld_snps.py \
    --diseases Alzheimers \
    --gwas_cluster_file_in  /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir/file_based_gwas//Alzheimers/gwas_clusters/gwas_cluster_with_10_snps_around_rs11055612.pickle \
    --gwas_cluster_file_out /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir/file_based_gwas//Alzheimers/gwas_clusters/gwas_cluster_with_10_snps_around_rs11055612.pvalues_from_file.pickle

	'''
	
	import postgap.Globals
	postgap.Globals.DATABASES_DIR = '/nfs/nobackup/ensembl/mnuhn/postgap/databases/'

	parser = argparse.ArgumentParser(description='Load pvalues for ld snps')
	parser.add_argument('--diseases', nargs='*', default=[])
	#parser.add_argument('--iris',     nargs='*', default=[])
	parser.add_argument('--gwas_cluster_file_in')
	parser.add_argument('--gwas_cluster_file_out')
	
	options = parser.parse_args()

	logger.info( "gwas_cluster_file_in   = " + options.gwas_cluster_file_in  )
	logger.info( "gwas_cluster_file_out  = " + options.gwas_cluster_file_out )
	#logger.info( "iris                   = " + pformat(options.iris) )
	logger.info( "diseases               = " + pformat(options.diseases)     )
	
	import pickle
	
	pickle_fh    = open(options.gwas_cluster_file_in, 'rb')
	gwas_cluster = pickle.load(pickle_fh)
	
	from postgap.Summarisers import summarise
	
	logger.info( "The Gwas cluster without substitutions from the file:\n" + summarise(gwas_cluster)     )

	diseases = options.diseases
	#iris     = options.iris
	iris     = []
	
	from postgap.FinemapIntegration.GwasIntegration import load_gwas_pvalues_from_file
	gwas_clusters_with_values_from_file, pvalues_were_loaded = load_gwas_pvalues_from_file([ gwas_cluster ], diseases, iris)

	# One cluster goes into load_gwas_pvalues_from_file, so only one should 
	# come out.
	#
	if len(gwas_clusters_with_values_from_file) != 1:
		raise Exception
	
	gwas_cluster_with_values_from_file = gwas_clusters_with_values_from_file[0]
	
	if pvalues_were_loaded:
		logging.info("The Gwas Cluster with substitutions from the file:\n" + summarise(gwas_cluster_with_values_from_file))
	else:
		logging.info("Gwas Cluster remained unchanged.")

	import os
	gwas_cluster_file_out = options.gwas_cluster_file_out
	output_directory = os.path.dirname(gwas_cluster_file_out)
	
	if not os.path.exists(output_directory):
		os.makedirs(output_directory)

	pickle_fh = open(gwas_cluster_file_out, 'w')
	pickle.dump(gwas_cluster_with_values_from_file, pickle_fh)
	pickle_fh.close
	
	logging.info("The gwas cluster has been written to " + gwas_cluster_file_out)
	logging.info("All done.");

if __name__ == "__main__":
	main()

