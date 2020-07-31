# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a script file to record memory consumption.
"""



import argparse
from argparse import RawTextHelpFormatter
import collections
import cPickle as pickle
import itertools as it
import json
import logging
import logging.config
import math
import numpy
import operator
import os
import os.path
import pprint
from pprint import pformat
import random
from subprocess import Popen, PIPE
import sys
import re
import tempfile

import postgap
import postgap.Cisreg
import postgap.EFO
import postgap.Ensembl_lookup
import postgap.Finemap
import postgap.FinemapIntegration
import postgap.Globals
import postgap.GWAS
import postgap.Integration
import postgap.LD
import postgap.Reg
import postgap.RegionFilter
import postgap.REST
from postgap.DataModel import *
from postgap.Utils import *

from guppy import hpy
h = hpy()


r2_cache = collections.defaultdict(dict)
r2_sets = set()
GRCh38_snp_locations= dict()
known_chroms = map(str, range(1,23)) + ['X','Y']

postgap.Globals.ALL_TISSUES = postgap.Integration.get_all_tissues()
#postgap.Globals.GWAS_SUMMARY_STATS_FILE = 'CAD_UKBIOBANK-chr22.tsv'
postgap.Globals.CLUSTER_FILE = '/homes/yalan/yalan/NNRtest/CAD_UKBIOBANK_200728124133/CAD_UKBIOBANK-chr9_clusters/GWAS_Cluster_9:17458027-17458027'
postgap.Globals.GTEx_path = '/nfs/production/panda/ensembl/funcgen/eqtl/GTEx.V6.88_38.cis.eqtls.h5'
postgap.Globals.SQLite_connection = '/nfs/production/panda/ensembl/funcgen/eqtl/GTEx.V6.88_38.cis.eqtls.h5.sqlite3'
postgap.Globals.SPECIES = 'Human'
postgap.Globals.DATABASES_DIR = '/homes/yalan/yalan/databases/'
postgap.Globals.PERFORM_BAYESIAN = True
postgap.Globals.KSTART = 1
postgap.Globals.KMAX = 2

debug = True
if debug: logging.getLogger().setLevel(logging.DEBUG)

population = 'EUR'
tissues = None
tissue_weights = None

efos = None
diseases = []
efo_iris = []
efo_iris = filter(lambda X: X is not None, (postgap.EFO.suggest(disease) for disease in diseases))
expanded_efo_iris = efo_iris

scipy = None
coords = None
rsID = None



# FinemapIntegration.py
def finemap_gwas_cluster(cluster, population):
	if len(cluster.ld_snps) == len(cluster.gwas_snps):
		ld_snps, ld_matrix, z_scores, betas = compute_ld_matrix(cluster, population)
	elif postgap.Globals.GWAS_SUMMARY_STATS_FILE is not None:
		ld_snps, ld_matrix, z_scores, betas = extract_z_scores_from_file(cluster, population)
	else:
		ld_snps, ld_matrix, z_scores, betas = impute_z_scores(cluster, population)
	## Define experiment label (serves for debugging logs)
	chrom = ld_snps[0].chrom
	start = min(ld_snp.pos for ld_snp in ld_snps)
	end = max(ld_snp.pos for ld_snp in ld_snps)
	sample_label = 'GWAS_Cluster_%s:%i-%i' % (chrom, start, end)
	## Define sample size: mean of max for each SNP
	sample_sizes = map(lambda gwas_snp: max(gwas_association.sample_size for gwas_association in gwas_snp.evidence), cluster.gwas_snps)
	sample_size = sum(sample_sizes) / len(sample_sizes)
	ld_snp_ids = [ld_snp.rsID for ld_snp in ld_snps]
	print('\nin finemap_gwas_cluster, will finemap:')
	print h.heap()
	## Compute posteriors
	configuration_posteriors = postgap.Finemap.finemap(
		z_scores	 = numpy.array(z_scores),
		beta_scores	= numpy.array(betas),
		cov_matrix	 = ld_matrix,
		n			= sample_size,
		labels		 = ld_snp_ids,
		sample_label = sample_label,
		kstart		 = postgap.Globals.KSTART,
		kmax		 = postgap.Globals.KMAX
	)
	print('\nin finemap_gwas_cluster, after finemap:')
	print h.heap()
	assert len(ld_snps) ==	ld_matrix.shape[0]
	assert len(ld_snps) ==	ld_matrix.shape[1]
	return GWAS_Cluster(cluster.gwas_snps, ld_snps, ld_matrix, z_scores, configuration_posteriors) 

def compute_ld_matrix(cluster, population):
	## Compute ld_matrix
	ld_snp_ids, ld_matrix = postgap.LD.get_pairwise_ld(cluster.ld_snps, population)
	## Update list of LD SNPs
	ld_snp_hash = dict((ld_snp.rsID, ld_snp) for index, ld_snp in enumerate(cluster.ld_snps))
	ld_snps = [ld_snp_hash[rsID] for rsID in ld_snp_ids]
	# Aggregate z_scores into a single vector
	z_scores = [0] * len(ld_snps)
	betas = [0] * len(ld_snps)
	gwas_snp_hash = dict((gwas_snp.snp.rsID, gwas_snp) for gwas_snp in cluster.gwas_snps)
	for index, ld_snp in enumerate(ld_snps):
		z_scores[index] = gwas_snp_hash[ld_snp.rsID].z_score
		betas[index] = gwas_snp_hash[ld_snp.rsID].beta
	assert len(ld_snps) ==	ld_matrix.shape[0]
	assert len(ld_snps) ==	ld_matrix.shape[1]
	return ld_snps, ld_matrix, z_scores, betas

def extract_z_scores_from_file(cluster, population):
	ld_snp_hash = dict((ld_snp.rsID, ld_snp) for index, ld_snp in enumerate(cluster.ld_snps))
	## Search ld snps in the GWAS summary stats file, and drop ld snps that cannot be found in the GWAS summary stats file from ld_snps
	proper_gwas_cluster = postgap.GWAS.GWAS_File().create_gwas_cluster_with_pvalues_from_file(gwas_cluster=cluster, gwas_data_file=postgap.Globals.GWAS_SUMMARY_STATS_FILE)
	## Extract z_scores for all the ld snps that can be found in the GWAS summary stats file
	ld_snp_results = dict((ld_snp.snp, (ld_snp.pvalue, ld_snp.beta)) for index, ld_snp in enumerate(proper_gwas_cluster.ld_snps))
	# Select all the found ld snps from ld_snps
	found_ld_snps = [ld_snp_hash[rsID] for rsID in ld_snp_results]
	## Compute ld_matrix
	ld_snp_ids, ld_matrix = postgap.LD.get_pairwise_ld(found_ld_snps, population)
	## Update the list of ld_snps
	ld_snps = [ld_snp_hash[rsID] for rsID in ld_snp_ids]
	z_scores = [z_score_from_pvalue(ld_snp_results[rsID][0], ld_snp_results[rsID][1]) for rsID in ld_snp_ids]
	betas = [ld_snp_results[rsID][1] for rsID in ld_snp_ids]
	assert len(ld_snps) ==	ld_matrix.shape[0]
	assert len(ld_snps) ==	ld_matrix.shape[1]
	return ld_snps, ld_matrix, z_scores, betas

def impute_z_scores(cluster, population):
	## Compute ld_matrix
	ld_snp_ids, ld_matrix = postgap.LD.get_pairwise_ld(cluster.ld_snps, population)
	## Update list of LD SNPs
	ld_snp_hash = dict((ld_snp.rsID, ld_snp) for index, ld_snp in enumerate(cluster.ld_snps))
	ld_snps = [ld_snp_hash[rsID] for rsID in ld_snp_ids]
	## Determine which SNPs are missing values
	gwas_snp_hash = dict((gwas_snp.snp.rsID, gwas_snp) for gwas_snp in cluster.gwas_snps)
	missing_indices = numpy.array([index for index, ld_snp in enumerate(ld_snps) if ld_snp.rsID not in gwas_snp_hash]).astype(int)
	known_z_scores = numpy.array([gwas_snp_hash[ld_snp.rsID].z_score for ld_snp in ld_snps if ld_snp.rsID in gwas_snp_hash])
	known_betas = numpy.array([gwas_snp_hash[ld_snp.rsID].beta for ld_snp in ld_snps if ld_snp.rsID in gwas_snp_hash])
	# Generate LD matrix of known values
	ld_matrix_known = numpy.delete(ld_matrix, missing_indices, axis=1)
	ld_matrix_known = numpy.delete(ld_matrix_known, missing_indices, axis=0)
	# Generate LD matrix of known SNPs to missing SNPs
	ld_matrix_k2m = ld_matrix[missing_indices, :]
	ld_matrix_k2m = numpy.delete(ld_matrix_k2m, missing_indices, axis=1)
	# Imputation
	shrink_lambda=0.1 # shrinkage factor, magic number
	ld_matrix_known_shrink = shrink_lambda *	numpy.diag(numpy.ones(ld_matrix_known.shape[0])) + (1-shrink_lambda) * ld_matrix_known
	ld_matrix_k2m_shrink = (1-shrink_lambda) * ld_matrix_k2m
	z_shrink_imputed = numpy.dot(numpy.dot(ld_matrix_k2m_shrink,numpy.linalg.pinv(ld_matrix_known_shrink, 0.0001)) , known_z_scores)
	beta_shrink_imputed = numpy.dot(numpy.dot(ld_matrix_k2m_shrink,numpy.linalg.pinv(ld_matrix_known_shrink, 0.0001)) , known_betas)
	# Aggregate z_scores into a single vector
	z_scores = []
	betas = []
	for i in ld_snps:
		z_scores.append(0)
		betas.append(0)
	for index, z_score, beta in zip(missing_indices, z_shrink_imputed, beta_shrink_imputed):
		z_scores[index] = z_score
		betas[index] = beta
	for index, ld_snp in enumerate(ld_snps):
		if ld_snp.rsID in gwas_snp_hash:
			z_scores[index] = gwas_snp_hash[ld_snp.rsID].z_score
			betas[index] = gwas_snp_hash[ld_snp.rsID].beta
	assert len(ld_snps) ==	ld_matrix.shape[0]
	assert len(ld_snps) ==	ld_matrix.shape[1]
	return ld_snps, ld_matrix, z_scores, betas

def compute_joint_posterior(cluster, associations):
	assert len(cluster.ld_snps) == cluster.ld_matrix.shape[0], (len(cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
	assert len(cluster.ld_snps) == cluster.ld_matrix.shape[1], (len(cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
	gene_tissue_snp_eQTL_hash = organise_eQTL_data(associations)
	return dict((gene, compute_gene_joint_posterior(cluster, gene, gene_tissue_snp_eQTL_hash[gene])) for gene in gene_tissue_snp_eQTL_hash)

def compute_gene_joint_posterior(cluster, gene, tissue_snp_eQTL_hash):
	assert len(cluster.ld_snps) == cluster.ld_matrix.shape[0], (len(cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
	assert len(cluster.ld_snps) == cluster.ld_matrix.shape[1], (len(cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
	return dict((tissue, compute_gene_tissue_joint_posterior(cluster, gene, tissue, tissue_snp_eQTL_hash[tissue])) for tissue in tissue_snp_eQTL_hash)

def compute_gene_tissue_joint_posterior(cluster, tissue, gene, eQTL_snp_hash):
	assert len(cluster.ld_snps) == cluster.ld_matrix.shape[0], (len(cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
	assert len(cluster.ld_snps) == cluster.ld_matrix.shape[1], (len(cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
	## Determine which SNPs are missing values
	missing_indices = numpy.array([index for index, ld_snp in enumerate(cluster.ld_snps) if ld_snp.rsID not in eQTL_snp_hash]).astype(int)
	known_z_scores = numpy.array([eQTL_snp_hash[ld_snp.rsID][0] for ld_snp in cluster.ld_snps if ld_snp.rsID in eQTL_snp_hash])
	known_betas = numpy.array([eQTL_snp_hash[ld_snp.rsID][1] for ld_snp in cluster.ld_snps if ld_snp.rsID in eQTL_snp_hash])
	assert all(beta is not None for beta in known_betas)
	assert len(known_z_scores) > 0
	assert len(missing_indices) != len(cluster.ld_snps), (missing_indices, known_z_scores)
	assert len(cluster.ld_snps) == cluster.ld_matrix.shape[0], (len(cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
	assert len(cluster.ld_snps) == cluster.ld_matrix.shape[1], (len(cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
	if len(missing_indices) > 0:
		# Generate LD matrix of known values
		ld_matrix_known = numpy.delete(cluster.ld_matrix, missing_indices, axis=1)
		ld_matrix_known = numpy.delete(ld_matrix_known, missing_indices, axis=0)
		assert ld_matrix_known.size > 0, (missing_indices, cluster.ld_matrix, ld_matrix_known)
		# Generate LD matrix of known SNPs to missing SNPs
		ld_matrix_k2m = cluster.ld_matrix[missing_indices, :]
		ld_matrix_k2m = numpy.delete(ld_matrix_k2m, missing_indices, axis=1)
		# Imputation
		shrink_lambda=0.1 # shrinkage factor, magic number
		ld_matrix_known_shrink = shrink_lambda *	numpy.diag(numpy.ones(ld_matrix_known.shape[0])) + (1-shrink_lambda) * ld_matrix_known
		assert ld_matrix_known_shrink.size > 0, (missing_indices, ld_matrix, ld_matrix_known)
		ld_matrix_k2m_shrink = (1-shrink_lambda) * ld_matrix_k2m
		z_shrink_imputed = numpy.dot(numpy.dot(ld_matrix_k2m_shrink,numpy.linalg.pinv(ld_matrix_known_shrink, 0.0001)) , known_z_scores)
		beta_shrink_imputed = numpy.dot(numpy.dot(ld_matrix_k2m_shrink,numpy.linalg.pinv(ld_matrix_known_shrink, 0.0001)) , known_betas)
		# Aggregate z_scores into a single vector
		z_scores = []
		betas = []
		for i in cluster.ld_snps:
			z_scores.append(0)
			betas.append(0)
		for index, z_score, beta in zip(missing_indices, z_shrink_imputed, beta_shrink_imputed):
			z_scores[index] = z_score
			betas[index] = beta
		for index, ld_snp in enumerate(cluster.ld_snps):
			if ld_snp.rsID in eQTL_snp_hash:
				z_scores[index] = eQTL_snp_hash[ld_snp.rsID][0]
				betas[index] = eQTL_snp_hash[ld_snp.rsID][1]
	else:
		z_scores = known_z_scores
		betas = known_betas
	## Define experiment label (serves for debugging logs)
	chrom = cluster.ld_snps[0].chrom
	start = min(ld_snp.pos for ld_snp in cluster.ld_snps)
	end = max(ld_snp.pos for ld_snp in cluster.ld_snps)
	sample_label = 'eQTL_Cluster_%s:%i-%i_%s' % (chrom, start, end, gene)
	## Compute posteriors
	eQTL_configuration_posteriors = postgap.Finemap.finemap(
		z_scores	 = numpy.array(z_scores),
		beta_scores	= numpy.array(betas),
		cov_matrix	 = cluster.ld_matrix,
		n			= 500, #TODO extract eQTL sample sizes
		kstart		 = postgap.Globals.KSTART,
		kmax		 = postgap.Globals.KMAX,
		sample_label = sample_label,
		labels		 = [ld_snp.rsID for ld_snp in cluster.ld_snps]
	)
	## Joint posterior
	sum_posteriors, config_sample = eQTL_configuration_posteriors.joint_posterior(cluster.gwas_configuration_posteriors)
	# Organise information into a hash 
	res = dict((config, config_sample.posterior[config_sample.configurations[config]]) for config in config_sample.configurations)
	res['_CLUSTER'] = sum_posteriors
	return res

def organise_eQTL_data(associations):
	res = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(float)))
	for association in associations:
		for evidence in association.cisregulatory_evidence:
			if evidence.source == 'GTEx':
				res[association.gene][evidence.tissue][association.snp.rsID] = (evidence.z_score, evidence.beta)
	return res

def sign(number):
	if number > 0:
		return 1
	elif number < 0:
		return -1
	else:
		return 0

def z_score_from_pvalue(p_value, direction):
	global scipy
	if scipy is None:
		import scipy
		import scipy.stats
	return -scipy.stats.norm.ppf(p_value/2) * sign(direction)


# Finemap.py
OneDConfigurationSample_prototype = collections.namedtuple(
	'OneDConfigurationSample', 
	[
		'configurations',
		'posterior',
		'log_BF',
		'configuration_size',
		'log_prior',
		'labels',
		'sample_label'
	]
)

class TranslationIndexException(Exception):
	pass

class OneDConfigurationSample(OneDConfigurationSample_prototype):
	def normalise_posteriors(self):
		max_logpp_unscaled = numpy.max(self.posterior)
		# Check if max_logpp_unscaled is so large that it overflows the exponential function 
		if numpy.exp(-max_logpp_unscaled) > 0:
			sum_calib = numpy.sum(numpy.exp(self.posterior-max_logpp_unscaled))
			assert not numpy.isinf(sum_calib), 'calibration going to infinity'
			normalisation_factor = numpy.exp(-max_logpp_unscaled) / sum_calib
			assert not numpy.isinf(normalisation_factor), 'normalisation going to infinity'
			return OneDConfigurationSample(
					configurations = self.configurations,
					posterior = numpy.exp(self.posterior) * normalisation_factor, 
					log_BF = self.log_BF,
					configuration_size = self.configuration_size,
					log_prior = self.log_prior,
					labels = self.labels,
					sample_label = self.sample_label
				)
		else:
			# If exponentiation is overflowed, reduce the posterior of non-maximal points to 0 
			posterior_max_points = numpy.where(self.posterior == max_logpp_unscaled, 1, 0)
			return OneDConfigurationSample(
					configurations = self.configurations,
					posterior = posterior_max_points / numpy.sum(posterior_max_points), 
					log_BF = self.log_BF,
					configuration_size = self.configuration_size,
					log_prior = self.log_prior,
					labels = self.labels,
					sample_label = self.sample_label
				)
	def marginals(self, singleton_count):
		singletons = range(singleton_count)
		configurations = dict((i, i) for i in singletons)
		configuration_size=numpy.ones(len(singletons))
		# compute marginal inclusion probability per SNP
		marginal = numpy.zeros(len(singletons))
		for configuration in self.configurations:
			for index in configuration:
				marginal[index] += self.posterior[self.configurations[configuration]]
		return OneDConfigurationSample(
			configurations = configurations,
			posterior = marginal,
			log_BF = None,
			configuration_size = configuration_size,
			log_prior = None,
			sample_label = self.sample_label
		)
	def multiple_marginals(self, singleton_count, kmax):
		iter_dict = {}
		for nc in range(1,(kmax+1)):
			#print(nc)
			config_iter = [config for config_size in range(nc, nc+1) for config in it.combinations(range(singleton_count), config_size)]
			len_config_iter=len(config_iter)
			configurations = dict((config_iter[i], i) for i in range(len_config_iter))
			configuration_size=numpy.ones(len_config_iter)*nc
			marginal_iter = numpy.zeros(len_config_iter)
			count=0
			input_keys=self.configurations.keys()
			while count < len(self.configurations):
				configuration =	input_keys[count]
				if len(configuration) < nc:
					count += 1
				else:
					index_tupel =	[config for config in it.combinations(configuration,nc)]
					for index in index_tupel:
						marginal_iter[configurations[index]] += self.posterior[self.configurations[configuration]]
						count += 1
			iter_out = OneDConfigurationSample(
				configurations = configurations,
				posterior = marginal_iter,
				log_BF = None,
				configuration_size = configuration_size,
				log_prior = None
			)
			iter_dict[nc] = iter_out
		return(iter_dict)
	def create_translation_index(self, from_labels, to_labels):
		translate = dict()
		for index_in_gwas, current_label in enumerate(from_labels):
			try:
				translate[index_in_gwas] = to_labels.index(current_label)
			except ValueError:
				translate[index_in_gwas] = None
		return translate
	def translate_configurations(self, translation_index, configurations):
		return dict(
			(self.translate_configuration(translation_index, configuration), configurations[configuration])
			for configuration in configurations
		)
	def translate_configuration(self, translation_index, configurations):
		return tuple(translation_index[ configuration_component ] for configuration_component in configurations)
	def joint_posterior(self, sample2):
		sample1_labels = self.labels
		sample2_labels = sample2.labels
		translate_from_sample2_to_sample1 = self.create_translation_index(
			from_labels = sample2_labels,
			to_labels	 = sample1_labels
		)
		translated_sample2_configurations = self.translate_configurations(
			translation_index = translate_from_sample2_to_sample1,
			configurations	= sample2.configurations
		)
		keys1 = set(self.configurations.keys())
		translated_keys2 = set(translated_sample2_configurations.keys())
		intersection = list(keys1 & translated_keys2)
		configurations = dict((configuration, index) for index, configuration in enumerate(intersection))
		ids1 = [self.configurations[configuration] for configuration in intersection]
		ids2 = [translated_sample2_configurations[configuration] for configuration in intersection]
		posterior1 = numpy.take(self.posterior, ids1)
		posterior2 = numpy.take(sample2.posterior, ids2)
		configuration_size = numpy.take(self.configuration_size, ids1)
		posterior = posterior1 * posterior2
		coloc_evidence = numpy.sum(posterior)
		return coloc_evidence, TwoDConfigurationSample(
				configurations = configurations,
				posterior = posterior,
				configuration_size = configuration_size,
				posterior1 = posterior1,
				posterior2 = posterior2,
				labels = self.labels,
				coloc_evidence = coloc_evidence,
				sample_1_label = self.sample_label,
				sample_2_label = sample2.sample_label,
		)
	def __str__(self):
		return self.string()
	def string(self, max_show_at_top = 5, max_show_at_end = 2):
		indentation = "	"
		total_configurations = len(self.configurations)
		if total_configurations >	max_show_at_top + max_show_at_end:
			num_show_at_top = max_show_at_top
			num_show_at_end = max_show_at_end
		else:
			num_show_at_top = total_configurations
			num_show_at_end = 0
		summary_lines = [ "A total of %i configurations have been investigated in %s:" % (total_configurations, self.sample_label) ]
		sorted_configurations = sorted(
			self.configurations, 
			key = lambda x: self.log_BF[self.configurations[x]],
			reverse = True
		)
		for configuration in sorted_configurations[:num_show_at_top]:
			summary_lines.append(indentation + self.configuration_string(configuration))
		if num_show_at_end > 0:
			summary_lines.append(indentation +	"...")
		for configuration in sorted_configurations[-num_show_at_end:]:
			summary_lines.append(indentation + self.configuration_string(configuration))
		return "\n".join(summary_lines)
	def configuration_string(self, configuration):
		index_of_configuration = configurations[configuration]
		posterior = self.posterior[index_of_configuration]
		prior	 = math.exp(self.log_prior[index_of_configuration])
		BF	= math.exp(log_BF[index_of_configuration])
		snps		= ', '.join([self.labels[position] for position in configuration])
		return "- The snp configuration ({}) has a prior probability of {:1.0%}. The posterior probability is {:.2e}. The Bayes factor is: {:.2e}".format(
				snps, prior, posterior, BF
			)


TwoDConfigurationSample_prototype = collections.namedtuple(
	'TwoDConfigurationSample_prototype', 
	[
		'configurations', 
		'posterior', 
		'configuration_size', 
		'posterior1', 
		'posterior2', 
		'labels', 
		'coloc_evidence',
		'sample_1_label',
		'sample_2_label'
	]
)

class TwoDConfigurationSample(TwoDConfigurationSample_prototype):
	def __str__(self):
		return self.string()
	def string(self, max_show_at_top = 5, max_show_at_end = 2):
		indentation = "	"
		total_configurations = len(self.configurations)
		if total_configurations >	max_show_at_top + max_show_at_end:
			num_show_at_top = max_show_at_top
			num_show_at_end = max_show_at_end
		else:
			num_show_at_top = total_configurations
			num_show_at_end = 0
		summary_lines = [ "A total of %i configurations have been investigated in %s and %s:" % (total_configurations, self.sample_1_label, self.sample_2_label) ]
		sorted_configurations = sorted(
			self.configurations, 
			key = lambda x: self.posterior[self.configurations[x]],
			reverse = True
		)
		for configuration in sorted_configurations[:num_show_at_top]:
			summary_lines.append(indentation + self.configuration_string(configuration))
		if num_show_at_end > 0:
			summary_lines.append(indentation +	"...")
		for configuration in sorted_configurations[-num_show_at_end:]:
			summary_lines.append(indentation + self.configuration_string(configuration))
		summary_lines.append("The coloc evidence is: {:2e}".format(self.coloc_evidence))
		return '\n'.join(summary_lines)
	def configuration_string(self, configuration):
		snp_labels = ", ".join([self.labels[position] for position in configuration])
		index = self.configurations[configuration]
		posterior	= self.posterior[index]
		posterior1 = self.posterior1[index]
		posterior2 = self.posterior2[index]
		return "The snp configuration ({}) had posterior probabilities of {:.2e} in {} and {:.2e} in {}. The joint posterior probability is {:.2e}.".format(
			snp_labels,
			float(posterior1),
			self.sample_1_label,
			float(posterior2),
			self.sample_2_label,
			float(posterior))

@profile
def finemap(z_scores, beta_scores, cov_matrix, n, labels, sample_label, kstart=1, kmax=5, corr_thresh=0.9, max_iter=100000, output="configuration", prior="independence_robust", v_scale=0.0025, g="BRIC", eigen_thresh=0.1, verbose=False):
	# Test inputs
	assert len(z_scores) == cov_matrix.shape[0], 'Covariance matrix has %i rows, %i expcted' % (cov_matrix.shape[0], len(z_scores))
	assert len(z_scores) == cov_matrix.shape[1], 'Covariance matrix has %i columns, %i expcted' % (cov_matrix.shape[0], len(z_scores))
	assert not kstart > kmax, 'Incorrect number of causal variants specified, kmax (%i) must be greater than kstart (%i)' % (kmax, kstart)
	assert not numpy.any(numpy.isnan(z_scores)), 'Missing values detected in z-scores'
	assert not numpy.any(numpy.isnan(cov_matrix)), 'Missing values detected in covariance matrix'
	cor_scores = beta_scores * numpy.sqrt(2*0.5*(0.5))
	print('\nin finemap_gwas_cluster -> finemap, will initialise compare_neighborhood:')
	print h.heap()
	# Initialise
	score_cache = dict()
	neighbourhood_cache = dict()
	configurations = [config for config_size in range(1, kstart + 1) for config in it.combinations(range(len(z_scores)), config_size)]
	results = compare_neighborhood(
		configs	= configurations, 
		z_scores = z_scores, 
		cor_scores = cor_scores,
		cov_matrix = cov_matrix, 
		kmax = kmax, 
		n = n, 
		score_cache = score_cache, 
		prior = prior, 
		v_scale = v_scale, 
		g = g, 
		labels=labels,
		sample_label = sample_label
	)
	# Simple search
	if kstart == kmax:
		res_out=results.normalise_posteriors()
		if output == "configuration":
			return res_out
		elif output == "marginal":
			return res_out.marginals(len(z_scores))
		else:
			assert False, "%s unkown" % {output}
	print('\nin finemap_gwas_cluster -> finemap, will SSS:')
	print h.heap()
	# shotgun stochastic search
	if kstart < kmax:
		p = results.normalise_posteriors().posterior
		current_config = configurations[numpy.random.choice(len(p), size=1, p=p)[0]]
		count = 1
		result_list = [results]
		while count < max_iter:
			# Generate new configs
			new_configs = create_neighborhood(current_config, len(z_scores), kstart, kmax, neighbourhood_cache)
			# Evaluate probabilities of these configs
			results_nh = compare_neighborhood(new_configs, z_scores, cor_scores, cov_matrix, kmax, n, score_cache, prior, corr_thresh,	v_scale=v_scale, g=g)
			# Add new entries into the results list
			result_list.append(results_nh)
			# Choose seed for next round among new configs
			prob = results_nh.normalise_posteriors().posterior
			current_config = new_configs[numpy.random.choice(len(new_configs), size=1, p=prob)[0]]
			# Chatter to stdout
			if verbose == True:
				print(str(numpy.round(float(count) / max_iter * 100)) + '% of the search done')
				print(current_config)
			# Keep count of sampled configs
			count += 1
		res_out = merge_samples(result_list, labels, sample_label).normalise_posteriors()
		print('\nin finemap_gwas_cluster -> finemap, after SSS:')
		print h.heap()
		if output == "configuration":
			return res_out
		elif output == "marginal":
			return res_out.marginals(len(z_scores))
		else:
			assert False, "%s unkown" % {output}

def create_neighborhood(current_config, m, kstart, kmax, neighbourhood_cache):
	if tuple(current_config) in neighbourhood_cache:
		return neighbourhood_cache[tuple(current_config)]	 
	current_size = len(current_config)
	new_configs = []
	if current_size <= kstart:
		# add move (m-k new possible config)
		for i in range(m):
			new_possible_config = numpy.append(current_config, i)
			new_possible_config.sort()
			new_possible_config = numpy.unique(new_possible_config)
			new_configs.append(new_possible_config)
		# swap move k(m-k)
		for i in range(current_size):
			for j in range(m):
				new_possible_config = numpy.delete(current_config, i)
				new_possible_config = numpy.append(new_possible_config, j)
				new_possible_config.sort()
				new_possible_config = numpy.unique(new_possible_config)
				new_configs.append(new_possible_config)
	if current_size > kstart and current_size < kmax:
		# delete move (k new possible config)
		for i in range(current_size):
			new_possible_config = numpy.delete(current_config, i)
			new_configs.append(new_possible_config)
		# add move (m-k new possible config)
		for i in range(m):
			new_possible_config = numpy.append(current_config, i)
			new_possible_config.sort()
			new_possible_config = numpy.unique(new_possible_config)
			new_configs.append(new_possible_config)
		# swap move k(m-k)
		for i in range(current_size):
			for j in range(m):
				new_possible_config = numpy.delete(current_config, i)
				new_possible_config = numpy.append(new_possible_config, j)
				new_possible_config.sort()
				new_possible_config = numpy.unique(new_possible_config)
				new_configs.append(new_possible_config)
	if current_size == kmax:
		# delete move (k new possible config)
		for i in range(current_size):
			new_possible_config = numpy.delete(current_config, i)
			new_configs.append(new_possible_config)
		# swap move k(m-k)
		for i in range(current_size):
			for j in range(m):
				new_possible_config = numpy.delete(current_config, i)
				new_possible_config = numpy.append(new_possible_config, j)
				new_possible_config.sort()
				new_possible_config = numpy.unique(new_possible_config)
				new_configs.append(new_possible_config)
	neighbourhood_cache[tuple(current_config)] = new_configs	
	return new_configs

def compare_neighborhood(configs, z_scores,	cor_scores, cov_matrix, kmax, n, score_cache, labels, sample_label, prior="independence_robust", corr_thresh=0.9, v_scale=0.0025, g="BRIC", eigen_thresh=0.1):
	nh_size = len(configs)
	configuration_size = numpy.array([len(configuration) for configuration in configs])
	log_prior = calc_logbinom(configuration_size, kmax, len(z_scores))
	log_BF = numpy.zeros(len(configs))
	i=0
	for configuration in configs:
		if tuple(configuration) in score_cache:
			log_BF[i] = score_cache[tuple(configuration)]
			i=i+1
			continue
		z_tuple = numpy.take(z_scores, configuration)
		cor_tuple = numpy.take(cor_scores, configuration)
		cov_tuple = cov_matrix[numpy.ix_(configuration, configuration)]
		if numpy.max(numpy.tril(cov_tuple,k=-1)) > corr_thresh:
			i = i+1
			continue			
		if numpy.min(numpy.tril(cov_tuple,k=-1)) < -corr_thresh:
			i = i+1
			continue
		tuple_size = len(z_tuple)
		if prior == "independence":
			v = numpy.eye(tuple_size) * v_scale
			log_BF[i] = calc_logBF(z_tuple, cov_tuple, v, n)
		elif prior == "independence_robust":
			v = numpy.eye(tuple_size) * v_scale
			log_BF[i] = calc_robust_logBF(z_tuple, cor_tuple, cov_tuple, v, n, eigen_thresh)
		elif prior == "gprior":
			log_BF[i] = calc_loggBF(z_tuple, cov_tuple, n, g=g)
		else:
			assert False, "%s is not one of independence or gprior" % (prior)
		score_cache[tuple(configuration)] = log_BF[i]
		i = i+1
	return OneDConfigurationSample(
			configurations = dict((tuple(config), index) for index, config in enumerate(configs)),
			posterior = (log_BF + log_prior),
			log_BF = log_BF,
			configuration_size = configuration_size,
			log_prior = log_prior,
			labels = labels,
			sample_label = sample_label
		)

def calc_logBF(z, cov, v, n):
	z = numpy.matrix(z)
	v = numpy.matrix(v)
	m = z.shape[1]
	try:
		coeff = 1. / math.sqrt(numpy.linalg.det((numpy.matrix(numpy.eye(m)) + n * numpy.matrix(v) * numpy.matrix(cov))))
	except ValueError:
		print cov
		eigen_vec=numpy.linalg.eig(numpy.matrix(cov))[0]
		Q=numpy.linalg.eig(numpy.matrix(cov))[1]
		eigen_pos = eigen_vec
		eigen_pos[eigen_vec<0] = 0
		cov=Q*numpy.diag(eigen_pos)*Q.transpose()
		coeff = 1. / math.sqrt(numpy.linalg.det((numpy.matrix(numpy.eye(m)) + n * numpy.matrix(v) * numpy.matrix(cov))))
	exponent = 0.5 * z * \
		numpy.matrix(numpy.linalg.pinv(((n * v).I + cov), 0.0001)) * z.T
	return numpy.array(numpy.log(coeff) + exponent)[0][0]

def calc_robust_logBF(z, cor, cov, v, n, eigen_thresh):
	z = numpy.matrix(z)
	cor = numpy.matrix(cor)
	v = numpy.matrix(v)
	m = z.shape[1]
	try: 
		coeff = 1. / math.sqrt(numpy.linalg.det((numpy.matrix(numpy.eye(m)) + n * numpy.matrix(v) * numpy.matrix(cov))))
	except ValueError:
		print cov
		eigen_vec=numpy.linalg.eig(numpy.matrix(cov))[0]
		Q=numpy.linalg.eig(numpy.matrix(cov))[1]
		eigen_pos = eigen_vec
		eigen_pos[eigen_vec<0] = 0
		cov=Q*numpy.diag(eigen_pos)*Q.transpose()
		coeff = 1. / math.sqrt(numpy.linalg.det((numpy.matrix(numpy.eye(m)) + n * numpy.matrix(v) * numpy.matrix(cov))))
	exponent = 0.5 * z * \
		numpy.matrix(numpy.linalg.pinv(((n * v).I + cov), 0.0001)) * z.T
	log_BF=numpy.array(numpy.log(coeff) + exponent)[0][0]
	if(m>1):
		if(check_eigenvals(cor,cov,eigen_thresh) == False): 
			log_BF=0.
	else: 
		log_BF=log_BF
	return log_BF

def calc_loggBF(z, cov, n, g="BRIC"):
	if(g=="BRIC"):
		 gp=numpy.max((len(z),n))
	if(g=="BIC"):
		 gp=n
	if(g=="RIC"):
		 gp=len(z)
	z = numpy.matrix(z)
	cov = numpy.matrix(cov)
	m = numpy.float(len(z))
	gp = numpy.float(gp)
	pinv = numpy.matrix(numpy.linalg.pinv(cov, 0.0001))
	coeff = (1 + gp)**(-m / 2)
	exponent = 0.5 * numpy.divide(gp, (gp + 1)) * z * pinv * z.T
	return numpy.array((math.log(coeff) + exponent))[0][0]

def calc_logbinom(subset_size, k, m):
	if k == 1:
		return numpy.zeros(m)
	else:
		p = float(1) / m
		p_binom = p**subset_size * (1 - p)**(m - subset_size)
		p_k = numpy.zeros(k - 1)
		for i in range(1, k):
			p_k[i - 1] = p**i * (1 - p)**(m - i)
		p_rescale = numpy.sum(p_k)
		p_out = p_binom / p_rescale
		return numpy.log(p_out)

def merge_samples(samples, labels, sample_label):
	configurations_old = dict((configuration, (sample, sample.configurations[configuration])) for sample in samples for configuration in sample.configurations)
	configurations = dict((configuration, index) for index, configuration in enumerate(configurations_old.keys()))
	posterior = numpy.zeros(len(configurations))
	configuration_size = numpy.zeros(len(configurations))
	log_BF = numpy.zeros(len(configurations))
	log_prior = numpy.zeros(len(configurations))
	for configuration in configurations:
		sample, old_index = configurations_old[configuration]
		new_index = configurations[configuration]
		posterior[new_index] = sample.posterior[old_index]
		configuration_size[new_index] = sample.configuration_size[old_index]
		log_BF[new_index] = sample.log_BF[old_index]
		log_prior[new_index] = sample.log_prior[old_index]
	return OneDConfigurationSample(
			configurations = configurations,
			posterior = posterior,
			log_BF = log_BF,
			configuration_size = configuration_size,
			log_prior = log_prior,
			labels = labels,
			sample_label = sample_label
		)
def check_eigenvals(cor,Sigma,eigen_thresh=0.1):
	assert cor.shape[1] == Sigma.shape[1], 'Covariance matrix has %i rows, %i expected' % (Sigma.shape[1], cor.shape[1])
	dim=cor.shape[1]
	S_tanh = numpy.matrix(numpy.zeros((dim+1, dim+1)))
	S_tanh[0,0] = 1
	S_tanh[1:dim+1,1:dim+1] = Sigma
	S_tanh[0,1:dim+1] = cor
	S_tanh[1:dim+1,0] = numpy.transpose(cor)
	return(numpy.min(numpy.linalg.eigvals(S_tanh)) > eigen_thresh)


# Integration.py
def cluster_to_genes(cluster, tissues, population):
	if postgap.Globals.PERFORM_BAYESIAN:
		assert len(cluster.ld_snps) == cluster.ld_matrix.shape[0], (len(cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
		assert len(cluster.ld_snps) == cluster.ld_matrix.shape[1], (len(cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
		# Obtain interaction data from LD snps
		associations = ld_snps_to_genes([snp for snp in cluster.ld_snps], tissues)
		gene_tissue_posteriors = compute_joint_posterior(cluster, associations)
		res = [GeneCluster_Association(gene=gene, score=None, collocation_posterior=gene_tissue_posteriors[gene], cluster=cluster, evidence=filter(lambda X: X.gene == gene, associations),	r2=None) for gene in gene_tissue_posteriors]
		logging.info("\tFound %i genes associated" % (len(res)))
	else:
		# Compute LD from top SNP
		top_gwas_hit = sorted(cluster.gwas_snps, key=lambda X: X.pvalue)[-1]
		ld = postgap.LD.get_lds_from_top_gwas(top_gwas_hit.snp, cluster.ld_snps, population=population)
		# Obtain interaction data from LD snps
		associations = ld_snps_to_genes([snp for snp in cluster.ld_snps if snp in ld], tissues)
		# Compute gene score
		gene_scores = dict(
		((association.gene, association.snp), (association, association.score * ld[association.snp]))
		for association in associations
		)
		if len(gene_scores) == 0: return []
		# OMIM exception
		max_score = max(X[1] for X in gene_scores.values())
		for gene, snp in gene_scores:
			if len(gene_to_phenotypes(gene)): gene_scores[(gene, snp)][1] = max_score
		# PICS score precomputed and normalised
		pics = PICS(ld, top_gwas_hit.pvalue)
		# Compute posterior
		res = [GeneCluster_Association(gene=gene, score=total_score(pics[snp], gene_scores[(gene, snp)][1]), collocation_posterior=None, cluster=cluster, evidence=gene_scores[(gene, snp)][:1], r2=ld[snp]) for (gene, snp) in gene_scores	if snp in pics]
		logging.info("\tFound %i genes associated around GWAS SNP %s" % (len(res), top_gwas_hit.snp.rsID))
	# Pick the association with the highest score
	return sorted(res, key=lambda X: X.score)

def ld_snps_to_genes(ld_snps, tissues):
	cisreg = cisregulatory_evidence(ld_snps, tissues) # Hash of hashes SNP => Gene => Cisregulatory_Evidence
	reg = regulatory_evidence(cisreg.keys(), tissues) # Hash: SNP => [ Regulatory_evidence ]
	return concatenate((create_SNP_GeneSNP_Associations(snp, reg[snp], cisreg[snp]) for snp in cisreg))

def cisregulatory_evidence(ld_snps, tissues):
	if postgap.Globals.Cisreg_adaptors == None:
		logging.info("Searching for cis-regulatory data on %i SNPs in all databases" % (len(ld_snps)))
		evidence = concatenate(source().run(ld_snps, tissues) for source in postgap.Cisreg.sources)
	else:
		logging.info("Searching for cis-regulatory data on %i SNPs in (%s)" % (len(ld_snps), ", ".join(postgap.Globals.Cisreg_adaptors)))
		evidence = concatenate(source().run(ld_snps, tissues) for source in postgap.Cisreg.get_filtered_subclasses(postgap.Globals.Cisreg_adaptors))
	flaten_evidence = []
	for evi in evidence:
		if type(evi)==list:
			for e in evi:
				flaten_evidence.append(e)
		else:
			flaten_evidence.append(evi)
	filtered_evidence = filter(lambda association: association is not None and association.gene is not None and association.gene.biotype == "protein_coding", flaten_evidence)
	# Group by snp, then gene:
	res = collections.defaultdict(lambda: collections.defaultdict(list))
	for association in filtered_evidence:
		assert type(association) is postgap.DataModel.Cisregulatory_Evidence, "association is Cisregulatory_Evidence"
		res[association.snp][association.gene].append(association)
	if postgap.Globals.Cisreg_adaptors == None:
		logging.debug(("Found %i cis-regulatory interactions in all databases" % (len(res))))
	else:
		logging.debug(("Found %i cis-regulatory interactions in (%s)" % (len(res), ", ".join(postgap.Globals.Cisreg_adaptors))))
	return res

def regulatory_evidence(snps, tissues):
	if postgap.Globals.Reg_adaptors == None:
		logging.info("Searching for regulatory data on %i SNPs in all databases" % (len(snps)))
		res = concatenate(source().run(snps, tissues) for source in postgap.Reg.sources)
		logging.info("Found %i regulatory SNPs among %i in all databases" % (len(res), len(snps)))
	else:
		logging.info("Searching for regulatory data on %i SNPs in (%s)" % (len(snps), ", ".join(postgap.Globals.Reg_adaptors)))
		res = concatenate(source().run(snps, tissues) for source in postgap.Reg.get_filtered_subclasses(postgap.Globals.Reg_adaptors))
		logging.info("Found %i regulatory SNPs among %i in (%s)" % (len(res), len(snps), ", ".join(postgap.Globals.Reg_adaptors)))
	# Group by SNP
	hash = collections.defaultdict(list)
	for hit in res:
		hash[hit.snp].append(hit)
	return hash

def create_SNP_GeneSNP_Associations(snp, reg, cisreg):
	intermediary_scores, gene_scores = compute_v2g_scores(reg, cisreg)
	rank = dict((score, index) for index, score in enumerate(sorted(gene_scores.values(), reverse=True)))
	return [GeneSNP_Association(gene=gene, snp=snp, cisregulatory_evidence=cisreg[gene], regulatory_evidence=reg, intermediary_scores=intermediary_scores[gene], score=gene_scores[gene], rank=rank[gene_scores[gene]] + 1) for gene in cisreg]

def compute_v2g_scores(reg, cisreg):
	intermediary_scores = dict()
	gene_scores = dict()
	for gene in cisreg:
		intermediary_scores[gene] = collections.defaultdict(int)
		seen = set()
		for evidence in cisreg[gene] + reg:
			if evidence.source not in seen or float(evidence.score) > intermediary_scores[gene][evidence.source]:
				intermediary_scores[gene][evidence.source] = float(evidence.score)
				seen.add(evidence.source)
			# VEP stats
			if evidence.source == 'VEP':
				intermediary_scores[gene]['VEP_count'] += 1
				intermediary_scores[gene]['VEP_sum'] += float(evidence.score) 
			if evidence.source == 'GTEx':
				intermediary_scores[gene][evidence.tissue] = float(evidence.score)
		# PCHiC
		intermediary_scores[gene]['PCHiC'] = min(intermediary_scores[gene]['PCHiC'], 1)
		# VEP
		if 'VEP' in intermediary_scores[gene]:
			intermediary_scores[gene]['VEP_mean'] = intermediary_scores[gene]['VEP_sum'] / intermediary_scores[gene]['VEP_count']
		# Weighted sum
		gene_scores[gene] = sum(intermediary_scores[gene][source] * postgap.Globals.EVIDENCE_WEIGHTS[source] for source in intermediary_scores[gene] if source in postgap.Globals.EVIDENCE_WEIGHTS)
	return intermediary_scores, gene_scores


# POSTGAP.py
def pretty_gene_output(associations):
	text = ['\t'.join(str(t) for t in [gene_association.gene.id, gene_association.cluster.gwas_configuration_posteriors.sample_label, snp_id, dict_posterior[(i,)], tissue, dict_posterior['_CLUSTER']]) for gene_association in associations for tissue, dict_posterior in gene_association.collocation_posterior.items() for i,snp_id in enumerate(gene_association.cluster.gwas_configuration_posteriors.labels)]
	return '\n'.join(text)

def pretty_output(associations, population):
	column_names = ['ld_snp_rsID', 'chrom', 'pos', 'GRCh38_chrom', 'GRCh38_pos', 'afr', 'amr', 'eas', 'eur', 'sas', 'gnomad', 'gnomad_sas', 'gnomad_oth', 'gnomad_asj', 'gnomad_nfe', 'gnomad_afr', 'gnomad_amr', 'gnomad_fin', 'gnomad_eas','gene_symbol', 'gene_id', 'gene_chrom', 'gene_tss', 'GRCh38_gene_chrom', 'GRCh38_gene_pos', 'disease_name', 'disease_efo_id', 'score', 'rank', 'r2', 'cluster_id', 'gwas_source', 'gwas_snp', 'gwas_pvalue', 'gwas_pvalue_description', 'gwas_odds_ratio', 'gwas_odds_ratio_ci_start', 'gwas_odds_ratio_ci_end', 'gwas_beta', 'gwas_size', 'gwas_pmid', 'gwas_study', 'gwas_reported_trait', 'ls_snp_is_gwas_snp', 'vep_terms', 'vep_sum', 'vep_mean'] + ["GTEx_" + tissue_name for tissue_name in postgap.Globals.ALL_TISSUES] + [source.display_name for source in postgap.Cisreg.sources + postgap.Reg.sources]
	if postgap.Globals.PERFORM_BAYESIAN: 
		column_names += [tissue_name + "_CLPP" for tissue_name in postgap.Globals.ALL_TISSUES]
	header = "\t".join(column_names).encode('utf-8')
	content = filter(lambda X: len(X) > 0, [pretty_cluster_association(association, population) for association in associations])
	return "\n".join([header] + content)

def pretty_cluster_association(association, population):
	results = genecluster_association_table(association, population)
	return "\n".join("\t".join([unicode(element).encode('utf-8') for element in row]) for row in results)

def genecluster_association_table(association, population):
	results = []
	snp_set = frozenset(association.cluster.ld_snps)
	if snp_set not in r2_sets:
		for snp1 in association.cluster.ld_snps:
			for snp2 in association.cluster.ld_snps:
				if snp1.rsID not in r2_cache or snp2.rsID not in r2_cache[snp1.rsID]:
					try: 
						ld_snp_ids, r_matrix = postgap.LD.get_pairwise_ld(association.cluster.ld_snps, population)
					except postgap.LD.UnitLDMatrixerror:
						break
					r2_sets.add(snp_set)
					r_index = dict((snp, index) for index, snp in enumerate(ld_snp_ids))
					for SNPA in ld_snp_ids:
						for SNPB in ld_snp_ids:
							r2_cache[SNPA][SNPB] = r_matrix.item((r_index[SNPA], r_index[SNPB]))**2
					break
	GRCh38_gene = postgap.Ensembl_lookup.get_ensembl_gene(association.gene.id, postgap.Ensembl_lookup.GRCH38_ENSEMBL_REST_SERVER)
	if GRCh38_gene is None:
		logging.info("%s not mapped onto GRCh38 - skipping" % association.gene.id)
		return []
	if GRCh38_gene.chrom not in known_chroms:
		logging.info("%s not on the principal GRCh38 assembly - skipping" % GRCh38_gene.chrom)
		return []
	GRCh38_snps = postgap.Ensembl_lookup.get_snp_locations([ld_snp.rsID for ld_snp in association.cluster.ld_snps if ld_snp.rsID not in GRCh38_snp_locations], postgap.Ensembl_lookup.GRCH38_ENSEMBL_REST_SERVER)
	for snp in GRCh38_snps:
		GRCh38_snp_locations[snp.rsID] = snp
	cluster_id = hash(json.dumps(association.cluster.gwas_snps))
	for gwas_snp in association.cluster.gwas_snps:
		for gwas_association in gwas_snp.evidence:
			pmid = clean_pmid(gwas_association.publication)
			if pmid is None:
				logging.info("PMID missing - skipping variant")
				continue
			for gene_snp_association in association.evidence:
				if gene_snp_association.snp.rsID not in GRCh38_snp_locations:
					logging.info("%s not on GRCh38 - skipping" % gene_snp_association.snp.rsID)
					continue
				if gene_snp_association.snp.chrom not in known_chroms or GRCh38_snp_locations[gene_snp_association.snp.rsID].chrom not in known_chroms:
					logging.info("%s not on main GRCh37 (%s) & GRCh38 assembly (%s) - skipping" % (gene_snp_association.snp.rsID, gene_snp_association.snp.chrom, GRCh38_snp_locations[gene_snp_association.snp.rsID].chrom))
					continue
				afr = 'N/A'
				amr = 'N/A'
				eas = 'N/A'
				eur = 'N/A'
				sas = 'N/A'
				gnomad = 'N/A'
				gnomad_sas = 'N/A'
				gnomad_oth = 'N/A'
				gnomad_asj = 'N/A'
				gnomad_nfe = 'N/A'
				gnomad_afr = 'N/A'
				gnomad_amr = 'N/A'
				gnomad_fin = 'N/A'
				gnomad_eas = 'N/A'
				vep_terms = []
				for evidence in gene_snp_association.cisregulatory_evidence:
					if evidence.source == "VEP":
						vep_terms += evidence.info['consequence_terms']
				if len(vep_terms) > 0:
					vep_terms = ",".join(list(set(vep_terms)))
				else:
					vep_terms = "N/A"
				for evidence in gene_snp_association.regulatory_evidence:
					if evidence.source == "VEP_reg":
						MAFs = evidence.info['MAFs']
						if MAFs is not None:
							if 'afr' in MAFs:
								afr = MAFs['afr']
							if 'amr' in MAFs:
								amr = MAFs['amr']
							if 'eas' in MAFs:
								eas = MAFs['eas']
							if 'eur' in MAFs:
								eur = MAFs['eur']
							if 'sas' in MAFs:
								sas = MAFs['sas']
							if 'gnomad' in MAFs:
								gnomad = MAFs['gnomad']
							if 'gnomad_sas' in MAFs:
								gnomad_sas = MAFs['gnomad_sas']
							if 'gnomad_oth' in MAFs:
								gnomad_oth = MAFs['gnomad_oth']
							if 'gnomad_asj' in MAFs:
								gnomad_asj = MAFs['gnomad_asj']
							if 'gnomad_nfe' in MAFs:
								gnomad_nfe = MAFs['gnomad_nfe']
							if 'gnomad_afr' in MAFs:
								gnomad_afr = MAFs['gnomad_afr']
							if 'gnomad_amr' in MAFs:
								gnomad_amr = MAFs['gnomad_amr']
							if 'gnomad_fin' in MAFs:
								gnomad_fin = MAFs['gnomad_fin']
							if 'gnomad_eas' in MAFs:
								gnomad_eas = MAFs['gnomad_eas']
							break
				if 'VEP_mean' in gene_snp_association.intermediary_scores:
					vep_mean = gene_snp_association.intermediary_scores['VEP_mean']
				else:
					vep_mean = 0
				if 'VEP_sum' in gene_snp_association.intermediary_scores:
					vep_sum = gene_snp_association.intermediary_scores['VEP_sum']
				else:
					vep_sum = 0
				tissue_eQTL_scores = []
				for tissue_name in postgap.Globals.ALL_TISSUES:
					if tissue_name in gene_snp_association.intermediary_scores:
						tissue_eQTL_scores.append(gene_snp_association.intermediary_scores[tissue_name])
					else:
						tissue_eQTL_scores.append(0)
				clpp = []
				if postgap.Globals.PERFORM_BAYESIAN:
					for tissue in postgap.Globals.ALL_TISSUES:
						if tissue in association.collocation_posterior:
							#clpp.append(sum(association.collocation_posterior[tissue][config] for config in association.collocation_posterior[tissue] if gene_snp_association.snp.rsID in config))
							clpp.append(0)
						else:
							clpp.append(0)
				r2_distance = read_pairwise_ld(gene_snp_association.snp, gwas_snp.snp)
				if r2_distance < 0.7:
					logging.info("%s LD from GWAS SNP < 0.7 - skipping" % gene_snp_association.snp.rsID)
					continue
				row = [
					gene_snp_association.snp.rsID, 
					gene_snp_association.snp.chrom, 
					gene_snp_association.snp.pos, 
					GRCh38_snp_locations[gene_snp_association.snp.rsID].chrom,
					GRCh38_snp_locations[gene_snp_association.snp.rsID].pos,
					afr,
					amr,
					eas,
					eur,
					sas,
					gnomad,
					gnomad_sas,
					gnomad_oth,
					gnomad_asj,
					gnomad_nfe,
					gnomad_afr,
					gnomad_amr,
					gnomad_fin,
					gnomad_eas,
					association.gene.name, 
					association.gene.id, 
					association.gene.chrom, 
					association.gene.tss, 
					GRCh38_gene.chrom,
					GRCh38_gene.tss,
					gwas_association.disease.name,
					re.sub(".*/", "", gwas_association.disease.efo),
					gene_snp_association.score, 
					gene_snp_association.rank,
					r2_distance,
					cluster_id,
					gwas_association.source,
					gwas_association.snp,
					gwas_association.pvalue,
					gwas_association.pvalue_description,
					gwas_association.odds_ratio,
					gwas_association.odds_ratio_ci_start,
					gwas_association.odds_ratio_ci_end,
					gwas_association.beta_coefficient,
					gwas_association.sample_size,
					pmid,
					gwas_association.study,
					gwas_association.reported_trait,
					int(gene_snp_association.snp.rsID == gwas_snp.snp.rsID),
					vep_terms,
					vep_sum,
					vep_mean
				]
				row += tissue_eQTL_scores
				row += [gene_snp_association.intermediary_scores[functional_source.display_name] for functional_source in postgap.Cisreg.sources + postgap.Reg.sources]
				if postgap.Globals.PERFORM_BAYESIAN:
					row += clpp
				results.append(row)
	return results

def clean_pmid(string):
	if re.match("PMID[0-9]+", string) is not None:
		return string
	elif re.match("[0-9]+", string) is not None:
		return "PMID" + string
	elif re.match("PUBMEDID:[0-9]+", string) is not None:
		return re.sub("PUBMEDID:", "PMID", string)
	else:
		return None

def read_pairwise_ld(snp1, snp2):
	if snp1.rsID == snp2.rsID:
		return 1
	if snp1.rsID in r2_cache and snp2.rsID in r2_cache:
		return r2_cache[snp1.rsID][snp2.rsID]
	else:
		return 0



print('\n\n\nafter loading all the dependencies:')
print h.heap()

logging.info("use cluster file, so skip previous steps and jump to cluster_to_genes (in gwas_snps_to_genes)")

if tissue_weights is None: tissue_weights = ['Whole_Blood']

cluster_fn = postgap.Globals.CLUSTER_FILE
infile = open(cluster_fn, 'rb')
cluster = pickle.load(infile)
infile.close()
print('\n\n\nafter loading cluster file in:')
print h.heap()

# Perform GWAS finemapping of the clusters
logging.info("\tperform GWAS finemapping for cluster %s" % (cluster_fn))

cluster = finemap_gwas_cluster(cluster, population)
print('\n\n\nafter GWAS finemapping:')
print h.heap()

res = cluster_to_genes(cluster, tissue_weights, population)
res = sorted(res, key=lambda X: X.score)
print('\n\n\nafter colocalisation:')
print h.heap()

logging.info("\tFound %i genes associated to cluster %s" % (len(res), cluster_fn))

pretty_gene_output(res)
print('\n\n\nafter organising the table in output2:')
print h.heap()

pretty_output(res, population)
print('\n\n\nafter organising the table in results:')
print h.heap()
