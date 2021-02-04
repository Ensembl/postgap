 # -*- coding: utf-8 -*-
"""
Spyder Editor

This is a script file to count SNPs in clusters
"""

import pickle

def get_n_snps_cluster(cluster):
	"""
	get 'n_ld_snps' and 'n_gwas_snps' from cluster file
	"""
	infile = open('CAD_UKBIOBANK-chr' + cluster.split(':')[0].split('_')[-1] + '_clusters/' + cluster, 'rb')
	clfile = pickle.load(infile)
	infile.close()
	return [len(clfile.ld_snps), len(clfile.gwas_snps)]
