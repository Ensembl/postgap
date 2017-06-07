#!/usr/bin/env python

# Run with this incantation:
#
# python -m unittest discover -s lib/t/ -p "*.py"

import unittest

class Blabladiblupp(unittest.TestCase):
 
	def setUp(self):
		pass
 
	def test_foooooo(self):
		
		import postgap.Globals
		postgap.Globals.DATABASES_DIR = '/nfs/nobackup/ensembl/mnuhn/postgap/databases/'
		pickle_file = 'lib/t/test-data/gwas_snps.merged_preclusters.pickle'
		
		pickle_fh = open(pickle_file, 'rb')
		
		import pickle
		clusters = pickle.load(pickle_fh)

		print type(clusters)
		print "Got %s clusters." % len(clusters)

		for x in range(len(clusters)):
			print "Cluster %s has %s members." % ( x, len(clusters[x].ld_snps) )

		import finemap.PostgapIntegation
		gwas_clusters_with_posteriors = finemap.PostgapIntegation.process_GWAS_Clusters(clusters)
		print "Got %i gwas clusters with posteriors." % len(gwas_clusters_with_posteriors)


if __name__ == '__main__':
	unittest.main()
