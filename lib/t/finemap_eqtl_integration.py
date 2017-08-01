#!/usr/bin/env python

# Run with this incantation:
#
# python -m unittest discover -s lib/t/ -p "*.py"

import unittest
import logging
import logging.config

class finemap_eqtl_integration(unittest.TestCase):
 
	def setUp(self):
		logging.config.fileConfig('configuration/logging.conf')
		pass
 
	def test_finemap_postgap_integration(self):
		
		logger = logging.getLogger(__name__)
		
		import postgap.Globals
		postgap.Globals.DATABASES_DIR = '/nfs/nobackup/ensembl/mnuhn/postgap/databases/'
		#pickle_file = 'SNP_GeneSNP_Associations.pickle'
		#pickle_file = 'SNP_GeneSNP_Associations.cisreg.pickle'
		pickle_file = 'output/EFO_0000203.pickle'
		
		pickle_fh = open(pickle_file, 'rb')
		
		import pickle
		snp_to_gene_and_evidence_dict = pickle.load(pickle_fh)

		logger.info(type(snp_to_gene_and_evidence_dict))
		logger.info("Got %s GeneSNP_Associations." % len(snp_to_gene_and_evidence_dict))
		
		import json
		import pprint
		
		pp = pprint.PrettyPrinter(indent=4, width=80)
		
		print "Hello world!"
		
		from postgap.DataModel import SNP, Gene, Cisregulatory_Evidence
		
		#
		# rehash snp_to_gene_and_evidence_dict which is a:
		#
		# dict[snp] -> dict( gene -> list (cisregulatory_evidence) )
		#
		# to tissue_gene_eqtl which is:
		#
		# dict[tissue][gene] -> list(cisregulatory_evidence)
		#
		
		import collections
		tissue_gene_eqtl = collections.defaultdict(generate_default)
		
		for snp in snp_to_gene_and_evidence_dict.keys():
			assert type(snp) is SNP, "snp is a SNP"
			
			gene_to_cisregulatory_evidence_list = snp_to_gene_and_evidence_dict[snp]
			
			for gene in gene_to_cisregulatory_evidence_list.keys():
				assert type(gene) is Gene, "gene is a Gene"
				
				cisregulatory_evidence_list = gene_to_cisregulatory_evidence_list[gene]
				assert type(cisregulatory_evidence_list) is list, "cisregulatory_evidence_list is a list"
				
				for cisregulatory_evidence in cisregulatory_evidence_list:					
					assert type(cisregulatory_evidence) is Cisregulatory_Evidence, "cisregulatory_evidence is Cisregulatory_Evidence"
					
					tissue_gene_eqtl[cisregulatory_evidence.tissue][gene].append(cisregulatory_evidence)
		#
		# Done creating tissue_gene_eqtl
		#
		
		# For each tissue and
		#    For each gene
		#       Run finemap
		#
		for tissue in tissue_gene_eqtl.keys():
			print tissue
			for gene in tissue_gene_eqtl[tissue].keys():
				print "    " + pp.pformat(gene)
				
				snps_correlated_to_current_tissue_and_gene = []
				
				for cisregulatory_evidence in tissue_gene_eqtl[tissue][gene]:
					snps_correlated_to_current_tissue_and_gene.append(cisregulatory_evidence.snp)
				
				# Finemap needs zscores.
				#
				# In order to compute zscores the betas are needed in addition 
				# to the pvalues.
				#
				# The betas are fetched from the rest server.
				#
				zscore_vector = []
				
				snp_hash = dict( (snp.rsID, snp) for snp in snps_correlated_to_current_tissue_and_gene)
				
				GTEx = postgap.Cisreg.GTEx()
				cisregulatory_evidence_with_betas = GTEx.gene_tissue_betas(gene, tissue, snp_hash)
				
				snp_to_beta = dict(
					(cisregulatory_evidence_with_beta.snp.rsID, cisregulatory_evidence_with_beta.beta) 
						for cisregulatory_evidence_with_beta in cisregulatory_evidence_with_betas
				)
				
				# The betas are in snp_to_beta now.
				
				print "snp_to_beta:"
				pp.pprint(snp_to_beta)
				
				# Compute vector of zscores.
				#
				for cisregulatory_evidence in tissue_gene_eqtl[tissue][gene]:
					pvalue = cisregulatory_evidence.pvalue
					zscore = compute_z_score_from_pvalue_and_sign(
						pvalue, 
						sign(snp_to_beta[cisregulatory_evidence.snp.rsID])
					)
					zscore_vector.append(zscore)
				
				print "zscore_vector:"
				pp.pprint(zscore_vector)
				
				print "snps_correlated_to_current_tissue_and_gene:"
				pp.pprint(snps_correlated_to_current_tissue_and_gene)
				
				(SNP_ids, r2_array) = postgap.LD.get_pairwise_ld(snps_correlated_to_current_tissue_and_gene)
				
				logger.info("SNPs:")
				logger.info("\n" + stringify_vector(SNP_ids))
				
				logger.info("zscores:")
				logger.info("\n" + stringify_vector(zscore_vector))
				
				logger.info("Matrix of pairwise linkage disequilibria:")
				logger.info("\n" + stringify_matrix(r2_array))
				
				logger.info("Running finemap")
				
				kstart = 1
				kmax   = 1
				max_iter = "Not used when kstart == kmax"

				import postgap.Finemap as sss
				finemap_posteriors = sss.finemap(
					z_scores   = zscore_vector,
					cov_matrix = r2_array,
					n          = len(r2_array),
					kstart     = kstart,
					kmax       = kmax,
					max_iter   = max_iter,
					prior      = "independence"
				)
				logger.info("Done running finemap")
				
				if finemap_posteriors is None:
					raise Exception ("ERROR: Didn't get posteriors!")


def stringify_vector(vector):
	import json
	return json.dumps(vector, indent=4)
	
def stringify_matrix(matrix):
	
	import numpy as np
	# https://docs.scipy.org/doc/numpy/reference/generated/numpy.set_printoptions.html
	np.set_printoptions(precision=3, suppress=True, linewidth=150)
	return str(np.matrix(matrix))



def generate_default():
	import collections
	return collections.defaultdict(list)

def sign(value):
	if value>0:
		return 1
	
	# None<0 is true, so must test for that separately.
	if value is not None and value<0:
		return -1
	
	if value==0:
		return 0

	raise Exception


def compute_z_score_from_pvalue_and_sign(pvalue, sign):

	from scipy.stats import norm
	
	if sign is None:
		return None

	# norm.ppf means:
	# Percent point function (inverse of cdf) at q of the given RV.
	#
	# See:
	# https://docs.scipy.org/doc/scipy-0.14.0/reference/stats.html
	
	z_score = - norm.ppf(pvalue/2) * sign
	
	return z_score

if __name__ == '__main__':
	unittest.main()
