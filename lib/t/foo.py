#!/usr/bin/env python

# Run with this incantation:
#
# python -m unittest discover -s lib/t/ -p "*.py"

import unittest
import logging
import logging.config

class Blabladiblupp(unittest.TestCase):
 
	def setUp(self):
		logging.config.fileConfig('configuration/logging.conf')
		pass
 
	def test_foooooo(self):
		
		logger = logging.getLogger(__name__)
		
		import postgap.Globals
		postgap.Globals.DATABASES_DIR = '/nfs/nobackup/ensembl/mnuhn/postgap/databases/'
		pickle_file = 'lib/t/test-data/gwas_snps.merged_preclusters.pickle'
		
		pickle_fh = open(pickle_file, 'rb')
		
		import pickle
		clusters = pickle.load(pickle_fh)

		logger.info(type(clusters))
		logger.info("Got %s clusters." % len(clusters))

		for x in range(len(clusters)):
			logger.info("Cluster %s has %s members." % ( x, len(clusters[x].ld_snps) ))

		import finemap.PostgapIntegation
		gwas_clusters_with_posteriors = finemap.PostgapIntegation.process_GWAS_Clusters(clusters)
		logger.info("Got %i gwas clusters with posteriors." % len(gwas_clusters_with_posteriors))


if __name__ == '__main__':
	unittest.main()
