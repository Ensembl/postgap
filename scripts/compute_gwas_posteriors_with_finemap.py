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

python scripts/compute_gwas_posteriors_with_finemap.py \
    --gwas_cluster_file      finemap/EFO_0000203/gwas_clusters/gwas_cluster_0.pickle \
    --output_posteriors_file finemap/EFO_0000203/gwas_clusters/gwas_cluster_0.posteriors.pickle

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
	
	pickle_fh       = open(options.gwas_cluster_file, 'rb')
	gwas_cluster    = pickle.load(pickle_fh)
	
	from pprint import pformat
	
	#for gwas_snp in gwas_cluster.gwas_snps:
		#for gwas_association in gwas_snp.evidence:
			##print json.dumps(gwas_association.rest_hash["_links"]["riskAlleles"]["href"])
			#print pformat(gwas_association)
			
	#sys.exit();
	
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
	
	gwas_posteriors = process_gwas_cluster(gwas_cluster, cluster_name = title)
	
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

	gwas_clusters = [ gwas_cluster ]

	logger.info(type(gwas_clusters))
	logger.info("Got %s gwas_clusters." % len(gwas_clusters))

	for x in range(len(gwas_clusters)):
		logger.info("Cluster %s has %s members." % ( x, len(gwas_clusters[x].ld_snps) ))

	from methods.GWAS_Cluster import compute_gwas_cluster_with_finemap_posteriors
	gwas_clusters_with_posteriors = compute_gwas_cluster_with_finemap_posteriors(gwas_cluster, cluster_name = cluster_name)
	
	if not(gwas_clusters_with_posteriors is None):
		logger.info("Got %i gwas clusters with posteriors." % len(gwas_clusters_with_posteriors))
	else:
		logger.warning("Got no gwas clusters with posteriors!")
	
	return gwas_clusters_with_posteriors.finemap_posteriors

if __name__ == "__main__":
	main()

