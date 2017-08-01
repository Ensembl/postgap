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
import logging
import logging.config
from pprint import pformat

logging.config.fileConfig('configuration/logging.conf')
logger = logging.getLogger(__name__)

def main():
	'''

python scripts/compute_gwas_posteriors_with_finemap.py \
    --gwas_cluster_file      finemap/EFO_0000203/gwas_clusters/gwas_cluster_0.pickle \
    --output_posteriors_file finemap/EFO_0000203/gwas_clusters/gwas_cluster_0.posteriors.pickle

	'''
	
	import postgap.Globals
	postgap.Globals.DATABASES_DIR = '/nfs/nobackup/ensembl/mnuhn/postgap/databases/'

	parser = argparse.ArgumentParser(description='Run finemap')
	parser.add_argument('--diseases', nargs='*', default=[])
	#parser.add_argument('--iris',     nargs='*', default=[])
	parser.add_argument('--gwas_cluster_file')
	parser.add_argument('--output_posteriors_file')
	
	options = parser.parse_args()
	
	logger.info( "gwas_cluster_file      = " + options.gwas_cluster_file      )
	logger.info( "output_posteriors_file = " + options.output_posteriors_file )
	#logger.info( "iris                   = " + pformat(options.iris) )
	logger.info( "diseases               = " + pformat(options.diseases) )
	
	#import sys
	#sys.exit(0)
	
	import pickle
	
	pickle_fh    = open(options.gwas_cluster_file, 'rb')
	gwas_cluster = pickle.load(pickle_fh)
	
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
	
	diseases = options.diseases
	#iris     = options.iris
	iris     = []
	
	from postgap.FinemapIntegration.GwasIntegration import load_gwas_pvalues_from_file
	gwas_clusters_with_values_from_file = load_gwas_pvalues_from_file([ gwas_cluster ], diseases, iris)
	
	# One cluster goes into load_gwas_pvalues_from_file, so only one should 
	# come out.
	#
	if len(gwas_clusters_with_values_from_file) != 1:
		raise Exception
	
	gwas_cluster_with_values_from_file = gwas_clusters_with_values_from_file[0]
	
	from postgap.FinemapIntegration.GWAS_Cluster import ZScoreComputationException
	from postgap.FinemapIntegration.GWAS_SNP     import snp_in_multiple_gwas_associations_exception
	
	try:
		
		gwas_posteriors = process_gwas_cluster(gwas_cluster_with_values_from_file, cluster_name = title)
		
	except ZScoreComputationException as z:
		
		#
		# "logger.info" no longer works here, no idea, why. "logging" does 
		# work somehow, so going with that.
		#
		logging.warning(str(z))
		logging.warning("Skipping this cluster.")
		return
	
	except snp_in_multiple_gwas_associations_exception as e:

		logging.warning(str(e))
		logging.warning("Skipping this cluster.")
		return
		
	logger.info( "Gwas posteriors have been computed:\n" + str(gwas_posteriors) )
	
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

def process_gwas_cluster(gwas_cluster, cluster_name="Unnamed cluster"):

	from postgap.FinemapIntegration.GWAS_Cluster import compute_gwas_cluster_with_finemap_posteriors
	gwas_clusters_with_posteriors = compute_gwas_cluster_with_finemap_posteriors(gwas_cluster, cluster_name = cluster_name)
	logger.info("Got %i gwas clusters with posteriors." % len(gwas_clusters_with_posteriors))
	return gwas_clusters_with_posteriors.finemap_posteriors

if __name__ == "__main__":
	main()

