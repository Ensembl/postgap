#! /usr/bin/env python

import logging
import logging.config
import sys
import json
import postgap.LD
import postgap.Globals
from postgap.DataModel import GWAS_Cluster
import pprint, pickle
import finemap.PostgapIntegation

def main():
	postgap.Globals.DATABASES_DIR = '/nfs/nobackup/ensembl/mnuhn/postgap/databases/'
	pickle_file = 'gwas_snps.merged_preclusters.pickle'
	
	pickle_fh = open(pickle_file, 'rb')
	clusters = pickle.load(pickle_fh)

	print type(clusters)
	print "Got %s clusters." % len(clusters)

	for x in range(len(clusters)):
		print "Cluster %s has %s members." % ( x, len(clusters[x].ld_snps) )

	gwas_clusters_with_posteriors = finemap.PostgapIntegation.process_GWAS_Clusters(clusters)
	print "Got %i gwas clusters with posteriors." % len(gwas_clusters_with_posteriors)

if __name__ == '__main__':
	main()

