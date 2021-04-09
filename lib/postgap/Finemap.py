#! /usr/bin/env python

"""

Copyright [1999-2018] EMBL-European Bioinformatics Institute

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

import scipy
import numpy
import math
import itertools as it
import collections
import logging

import postgap.Globals
from postgap.DataModel import *

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
	'''
			Stores finemapping data on a set of configurations:
			configurations: dict(tuple => int); assigns to each configuration (tuple of indices) an index into the vectors below
			posterior: numpy.array (1D)
			log_BF: numpy.array (1D)
			log_prior: numpy.array (1D)
			configuration_size: numpy.array (1D), number of SNPs in the corresponding vector
	'''

	def normalise_posteriors(self):
		'''
				Return identical OneDConfigurationSample but with updated non-log posteriors that add up to 1
				Arg1: OneDConfigurationSample
				Returntype: OneDConfigurationSample
		'''
		self.posterior[self.posterior > 700.] = 700.
		max_logpp_unscaled = numpy.max(self.posterior)
		assert not numpy.isinf(numpy.exp(max_logpp_unscaled)), 'max BF going to infinity'
		sum_calib = numpy.sum(numpy.exp(self.posterior-max_logpp_unscaled)) * numpy.exp(max_logpp_unscaled)
		assert not numpy.isinf(sum_calib), 'calibration going to infinity'
		return OneDConfigurationSample(
			configurations=self.configurations,
			posterior=numpy.exp(self.posterior) / sum_calib,
			log_BF=self.log_BF,
			configuration_size=self.configuration_size,
			log_prior=self.log_prior,
			labels=self.labels,
			sample_label=self.sample_label
		)

	def marginals(self, singleton_count):
		'''
				Return OneDConfigurationSample where the configuration posteriors are reduced to marginals
				on the positions
				Arg2: OneDConfigurationSample
				Returntype: OneDConfigurationSample
		'''
		singletons = range(singleton_count)
		configurations = dict((i, i) for i in singletons)
		configuration_size = numpy.ones(len(singletons))

		# compute marginal inclusion probability per SNP
		marginal = numpy.zeros(len(singletons))
		for configuration in self.configurations:
			for index in configuration:
				marginal[index] += self.posterior[self.configurations[configuration]]

		return OneDConfigurationSample(
			configurations=configurations,
			posterior=marginal,
			log_BF=None,
			configuration_size=configuration_size,
			log_prior=None,
			labels=self.labels,
			sample_label=self.sample_label
		)

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
		return tuple(translation_index[configuration_component] for configuration_component in configurations)

	def joint_posterior(self, sample2):
		'''
				Main function for fine-mapping of two traits
				Arg1: OneDConfigurationSample
				Arg2: OneDConfigurationSample
				Returntype: coloc_evidence, TwoDConfigurationSample
		'''

		sample1_labels = self.labels
		sample2_labels = sample2.labels

		translate_from_sample2_to_sample1 = self.create_translation_index(
			from_labels=sample2_labels,
			to_labels=sample1_labels
		)

		translated_sample2_configurations = self.translate_configurations(
			translation_index=translate_from_sample2_to_sample1,
			configurations=sample2.configurations
		)

		keys1 = set(self.configurations.keys())
		translated_keys2 = set(translated_sample2_configurations.keys())

		intersection = list(keys1 & translated_keys2)
		configurations = dict((configuration, index)
							  for index, configuration in enumerate(intersection))
		ids1 = [self.configurations[configuration]
				for configuration in intersection]
		ids2 = [translated_sample2_configurations[configuration]
				for configuration in intersection]

		posterior1 = numpy.take(self.posterior, ids1)
		posterior2 = numpy.take(sample2.posterior, ids2)
		configuration_size = numpy.take(self.configuration_size, ids1)
		posterior = posterior1 * posterior2
		coloc_evidence = numpy.sum(posterior)

		return coloc_evidence, TwoDConfigurationSample(
			configurations=configurations,
			posterior=posterior,
			configuration_size=configuration_size,
			posterior1=posterior1,
			posterior2=posterior2,
			labels=self.labels,
			coloc_evidence=coloc_evidence,
			sample_1_label=self.sample_label,
			sample_2_label=sample2.sample_label,
		)

	def __str__(self):
		return self.string()

	def string(self, max_show_at_top=5, max_show_at_end=2):
		indentation = "	"
		total_configurations = len(self.configurations)
		if total_configurations > max_show_at_top + max_show_at_end:
			num_show_at_top = max_show_at_top
			num_show_at_end = max_show_at_end
		else:
			num_show_at_top = total_configurations
			num_show_at_end = 0

		summary_lines = ["A total of %i configurations have been investigated in %s:" % (total_configurations, self.sample_label)]

		sorted_configurations = sorted(
			self.configurations,
			key=lambda x: self.log_BF[self.configurations[x]],
			reverse=True
		)

		summary_lines += [indentation + self.configuration_string(configuration) for configuration in sorted_configurations[:num_show_at_top]]

		if num_show_at_end > 0:
			summary_lines.append(indentation + "...")

		summary_lines += [indentation + self.configuration_string(configuration) for configuration in sorted_configurations[-num_show_at_end:]]

		return "\n".join(summary_lines)

	def configuration_string(self, configuration):
		index_of_configuration = self.configurations[configuration]
		posterior = self.posterior[index_of_configuration]
		prior = math.exp(self.log_prior[index_of_configuration])
		BF = math.exp(self.log_BF[index_of_configuration])
		snps = ', '.join([self.labels[position] for position in configuration])
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
'''
	Stores colocalisation data on a set of configurations:
	configurations: dict(tuple => int); assigns to each configuration (tuple of indices) an index into the vectors below
	posterior: numpy.array (1D), joint posterior
	configuration_size: numpy.array (1D), number of SNPs in the corresponding vector
	posterior1: numpy.array (1D)
	log_BF1: numpy.array (1D)
	log_prior1: numpy.array (1D)
	posterior2: numpy.array (1D)
	log_BF2: numpy.array (1D)
	log_prior2: numpy.array (1D)
'''


class TwoDConfigurationSample(TwoDConfigurationSample_prototype):
	def __str__(self):
		return self.string()

	def string(self, max_show_at_top=5, max_show_at_end=2):
		indentation = "	"
		total_configurations = len(self.configurations)
		if total_configurations > max_show_at_top + max_show_at_end:
			num_show_at_top = max_show_at_top
			num_show_at_end = max_show_at_end
		else:
			num_show_at_top = total_configurations
			num_show_at_end = 0

		summary_lines = ["A total of %i configurations have been investigated in %s and %s:" % (total_configurations, self.sample_1_label, self.sample_2_label)]

		sorted_configurations = sorted(
			self.configurations,
			key=lambda x: self.posterior[self.configurations[x]],
			reverse=True
		)

		summary_lines += [indentation + self.configuration_string(configuration) for configuration in sorted_configurations[:num_show_at_top]]

		if num_show_at_end > 0:
			summary_lines.append(indentation + "...")

		summary_lines += [indentation + self.configuration_string(configuration) for configuration in sorted_configurations[-num_show_at_end:]]

		summary_lines.append("The coloc evidence is: {:2e}".format(self.coloc_evidence))

		return '\n'.join(summary_lines)

	def configuration_string(self, configuration):
		snp_labels = ", ".join([self.labels[position] for position in configuration])
		index = self.configurations[configuration]
		return "The snp configuration ({}) had posterior probabilities of {:.2e} in {} and {:.2e} in {}. The joint posterior probability is {:.2e}.".format(
			snp_labels,
			float(self.posterior1[index]),
			self.sample_1_label,
			float(self.posterior2[index]),
			self.sample_2_label,
			float(self.posterior[index]))


def finemap_v1(z_scores, beta_scores, cov_matrix, n, labels, sample_label, lambdas, mafs, annotations, kstart=1, kmax=1, corr_thresh=0.9, max_iter=100, prior="independence_robust", v_scale=0.0025, g="BRIC", eigen_thresh=0.1, verbose=False, isGWAS=False):
	'''
			Main function for fine-mapping using stochastic search for one trait #
			Arg1: z_scores: numpy.array
			Arg2: beta_scores: numpy.array
			Arg3: cov_matrix numpy.array, correlation-structure (TwoD) array
			Arg4: n: int, sample size
			Arg5: list of strings, snp rsIDs
			Arg6: name of cluster (arbitrary) for easier reporting
			Arg7: numpy.array (1D), annotation effect size which is estimated by MLE
			Arg8: numpy.array (1D), minor allele frequency at cluster
			Arg9: numpy.array (2D), functional annotation matrix
			Arg kstart: int, full exploration of sets with #kstart causal variants
			Arg kmax: int, maximum number of causal variants
			Arg corr_thresh: int, excluding configurations with correlation > corr_thresh
			Arg max_iter: int, iterations of stochastic search
			Arg prior= string ("independence", "independence_robust" or "gprior") independence_robust does filter out configurations with possible mis-matches between the correlation matrix and effect size direction
			Arg v_scale = float, prior variance of the independence prior, recommended 0.05 **2 (following FINEMAP, Benner et al 2016)
			Arg g: string, g-parameter of the g-prior, recommended g="BRIC" g=max(n,#SNPs**2), other options g="BIC" where g=n (Bayes Information Criterion), or  g="RIC" where g=#SNPs**2 (Risk Inflation Criterion) see Mixtures of g Priors for Bayesian Variable Selection Liang et al 2008
			Arg verbose = False: boolean, print the progress of the stochastic search
			Returntype: OneDConfigurationSample
	'''
	logging.debug("Finemap v1")
	
	# Test inputs
	assert len(z_scores) == cov_matrix.shape[0], 'Covariance matrix has %i rows, %i expected' % (cov_matrix.shape[0], len(z_scores))
	assert len(z_scores) == cov_matrix.shape[1], 'Covariance matrix has %i columns, %i expected' % (cov_matrix.shape[0], len(z_scores))
	assert not kstart > kmax, 'Incorrect number of causal variants specified, kmax (%i) must be greater than kstart (%s)' % (kmax, kstart)
	assert not numpy.any(numpy.isnan(z_scores)), 'Missing values detected in z-scores'
	assert not numpy.any(numpy.isnan(cov_matrix)), 'Missing values detected in covariance matrix'

	# compute the correlation from the aligned beta
	# note this is an approximation with maf = 0.5
	# these correlations are NOT used for inference,
	# only for qc if the correlation between trait and SNP is aligned with the SNP x SNP correlation matrix
	#gwas_cor_aligned = gwas_b_v2*allele_flip * numpy.sqrt(2*gwas_maf_hg37_v2*(1-gwas_maf_hg37_v2))

	cor_scores = beta_scores * numpy.sqrt(2*0.5*(0.5))

	# Initialise
	score_cache = dict()
	neighbourhood_cache = dict()
	configurations = [config for config_size in range(1, kstart + 1) for config in it.combinations(range(len(z_scores)), config_size)]
	results = compare_neighborhood(
		configs=configurations,
		z_scores=z_scores,
		cor_scores=cor_scores,
		cov_matrix=cov_matrix,
		kmax=kmax,
		n=n,
		score_cache=score_cache,
		labels=labels,
		sample_label=sample_label,
		annotations=annotations,
		lambdas=lambdas,
		prior=prior,
		v_scale=v_scale,
		g=g,
	)

	# shotgun stochastic search
	if kstart < kmax:
		# Count the amount of all the possible configurations
		n_all_pos_config = sum(scipy.special.comb(len(z_scores), r, exact=True) for r in range(kstart, kmax+1))
		
		p = results.normalise_posteriors().posterior
		current_config = configurations[numpy.random.choice(len(p), size=1, p=p)[
			0]]
		count = 1
		while count < max_iter:
			# Generate new configs
			new_configs = create_neighborhood(current_config, len(
				z_scores), kstart, kmax, neighbourhood_cache)

			# Evaluate probabilities of these configs
			results_nh = compare_neighborhood(configs=new_configs,
							  z_scores=z_scores,
							  cor_scores=cor_scores,
							  cov_matrix=cov_matrix,
							  kmax=kmax,
							  n=n,
							  score_cache=score_cache,
							  labels=labels,
							  sample_label=sample_label,
							  annotations=annotations,
							  lambdas=lambdas,
							  prior=prior,
							  corr_thresh=corr_thresh,
							  v_scale=v_scale,
							  g=g)

			# Add new entries into the results object
			results = merge_samples(results, results_nh, labels, sample_label)

			# if all the possible configurations have been searched, stop iterations
			if len(results.configurations) == n_all_pos_config:
				logging.info(str(n_all_pos_config) + ' possible configurations have been searched within ' + str(count) + ' iterations; breaking.')
				break

			# Choose seed for next round among new configs
			prob = results_nh.normalise_posteriors().posterior
			current_config = new_configs[numpy.random.choice(len(new_configs), size=1, p=prob)[0]]

			# Chatter to stdout
			logging.debug(str(numpy.round(float(count) / max_iter * 100)) + '% of the search done')
			logging.debug(current_config)

			# Keep count of sampled configs
			count += 1

	return results.normalise_posteriors()


def finemap_v2(z_scores, beta_scores, cov_matrix, n, labels, sample_label, lambdas, mafs, annotations, kstart=1, kmax=1, corr_thresh=0.9, max_iter=100, prior="independence_robust", v_scale=0.0025, g="BRIC", eigen_thresh=0.1, verbose=False, isGWAS=False):
	'''
			Main function for fine-mapping using stochastic search for one trait #
			Arg1: z_scores: numpy.array
			Arg2: beta_scores: numpy.array
			Arg3: cov_matrix numpy.array, correlation-structure (TwoD) array
			Arg4: n: int, sample size
			Arg5: list of strings, snp rsIDs
			Arg6: name of cluster (arbitrary) for easier reporting
			Arg7: numpy.array (1D), annotation effect size which is estimated by MLE
			Arg8: numpy.array (1D), minor allele frequency at cluster
			Arg9: numpy.array (2D), functional annotation matrix
			Arg kstart: int, full exploration of sets with #kstart causal variants
			Arg kmax: int, maximum number of causal variants
			Arg corr_thresh: int, excluding configurations with correlation > corr_thresh
			Arg max_iter: int, iterations of stochastic search
			Arg prior= string ("independence", "independence_robust" or "gprior") independence_robust does filter out configurations with possible mis-matches between the correlation matrix and effect size direction
			Arg v_scale = float, prior variance of the independence prior, recommended 0.05 **2 (following FINEMAP, Benner et al 2016)
			Arg g: string, g-parameter of the g-prior, recommended g="BRIC" g=max(n,#SNPs**2), other options g="BIC" where g=n (Bayes Information Criterion), or  g="RIC" where g=#SNPs**2 (Risk Inflation Criterion) see Mixtures of g Priors for Bayesian Variable Selection Liang et al 2008
			Arg verbose = False: boolean, print the progress of the stochastic search
			Returntype: OneDConfigurationSample
	'''
	logging.debug("Finemap v2")

	# Test inputs
	assert len(z_scores) == cov_matrix.shape[0], 'Covariance matrix has %i rows, %i expected' % (cov_matrix.shape[0], len(z_scores))
	assert len(z_scores) == cov_matrix.shape[1], 'Covariance matrix has %i columns, %i expected' % (cov_matrix.shape[0], len(z_scores))
	assert not kstart > kmax, 'Incorrect number of causal variants specified, kmax (%i) must be greater than kstart (%s)' % (kmax, kstart)
	assert not numpy.any(numpy.isnan(z_scores)), 'Missing values detected in z-scores'
	assert not numpy.any(numpy.isnan(cov_matrix)), 'Missing values detected in covariance matrix'

	# compute the correlation from the aligned beta
	# note this is an approximation with maf = 0.5
	# these correlations are NOT used for inference,
	# only for qc if the correlation between trait and SNP is aligned with the SNP x SNP correlation matrix
	#gwas_cor_aligned = gwas_b_v2*allele_flip * numpy.sqrt(2*gwas_maf_hg37_v2*(1-gwas_maf_hg37_v2))

	cor_scores = beta_scores * numpy.sqrt(2*0.5*(0.5))

	# Initialise
	score_cache = dict()
	neighbourhood_cache = dict()
	configurations = [config for config_size in range(1, kstart + 1) for config in it.combinations(range(len(z_scores)), config_size)]
	results = compare_neighborhood(
		configs=configurations,
		z_scores=z_scores,
		cor_scores=cor_scores,
		cov_matrix=cov_matrix,
		kmax=kmax,
		n=n,
		score_cache=score_cache,
		labels=labels,
		sample_label=sample_label,
		annotations=annotations,
		lambdas=lambdas,
		prior=prior,
		v_scale=v_scale,
		g=g
	)

	# shotgun stochastic search
	if kstart < kmax:
		# Count the amount of all the possible configurations
		n_all_pos_config = sum(scipy.special.comb(len(z_scores), r, exact=True) for r in range(kstart, kmax+1))

		p = results.normalise_posteriors().posterior
		current_config = configurations[numpy.random.choice(len(p), size=1, p=p)[0]]
		count = 1
		while count < max_iter:
			# Generate new configs
			new_configs = create_neighborhood(current_config, len(z_scores), kstart, kmax, neighbourhood_cache)

			# EM Algorithm
			# E-step
			# compute PIPs within the list
			results_temp = compare_neighborhood(configs=results.configurations,
									z_scores=z_scores,
									cor_scores=cor_scores,
									cov_matrix=cov_matrix,
									kmax=kmax,
									n=n,
									score_cache=score_cache,
									labels=labels,
									sample_label=sample_label,
									annotations=annotations,
									lambdas=lambdas,
									prior=prior,
									corr_thresh=corr_thresh,
									v_scale=v_scale,
									g=g)
			PIPs = results_temp.normalise_posteriors().marginals(len(z_scores)).posterior

			# Q cost function
			def Q_function(lambdas_):
				'''
					compute Q cost function in E-step
					Arg1: numpy.array (1D), annotation effect size which will be estimated by EM
					Returntype: float, negative Q value for minimization
				'''
				Q = 0
				for i in range(len(z_scores)):
					annot = annotations[:, i]
					logistic_const = (1 + numpy.exp(-sum(annot * lambdas)))
					logistic_param = (1 + numpy.exp(-sum(annot * lambdas_)))
					Q = Q + ((
						PIPs[i]/(1.-PIPs[i]) * 1./logistic_const ) 
						* numpy.log(1./logistic_param) 
						+ (1. - 1./logistic_const) 
						* numpy.log(1. - 1./logistic_param))
				return -Q

			# M-step
			# Estimate parameters by maximising Q cost function
			bounds_list = [(-30, 30) for i in range(len(lambdas))]
			result_Q = scipy.optimize.minimize(
				# Q_function, lambdas, method='L-BFGS-B')
				Q_function, lambdas, method='L-BFGS-B', bounds=bounds_list)
			if result_Q.success:
				# lambdas = 1 / (1 + numpy.exp(-result_Q.x))
				lambdas = 1.0 / (1.0 + numpy.exp(-result_Q.x))

			# Evaluate probabilities of these configs
			results_nh = compare_neighborhood(configs=new_configs,
								  z_scores=z_scores,
								  cor_scores=cor_scores,
								  cov_matrix=cov_matrix,
								  kmax=kmax,
								  n=n,
								  score_cache=score_cache,
								  labels=labels,
								  sample_label=sample_label,
								  annotations=annotations,
								  lambdas=lambdas,
								  prior=prior,
								  corr_thresh=corr_thresh,
								  v_scale=v_scale,
								  g=g)

			# Add new entries into the results object
			results = merge_samples(results, results_nh, labels, sample_label)

			# if all the possible configurations have been searched, stop iterations
			if len(results.configurations) == n_all_pos_config:
				logging.info(str(n_all_pos_config) + ' possible configurations have been searched within ' + str(count) + ' iterations; breaking.')
				break

			# Choose seed for next round among new configs
			prob = results_nh.normalise_posteriors().posterior
			current_config = new_configs[numpy.random.choice(len(new_configs), size=1, p=prob)[0]]

			# Chatter to stdout
			logging.debug(str(numpy.round(float(count) / max_iter * 100)) +'% of the search done')
			logging.debug(current_config)

			# Keep count of sampled configs
			count += 1

	return results.normalise_posteriors()

def create_neighborhood(current_config, m, kstart, kmax, neighbourhood_cache):
	'''
			Defines a new set of possible moves starting from current config
			1. add
			2. swap
			3. delete
			Arg1: np.array (a possible starting causal configuration)
			Arg2: int, the size of the locus
			Arg3: int
			Arg4: int
			Returntype: array[arrays], other possible causal configurations
	'''
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


def compare_neighborhood(configs, z_scores,  cor_scores, cov_matrix, kmax, n, score_cache, labels, sample_label, annotations, lambdas, prior="independence_robust", corr_thresh=0.9, v_scale=0.0025, g="BRIC", eigen_thresh=0.1):
	'''
			Compare the moves with respect to the unscaled log posterior probability
			Arg1: array of arrays
			Arg2: numpy.array (OneD)
			Arg3: numpy.array (OneD)
			Arg4: numpy.array (TwoD)
			Arg5: kmax: int, maximum number of causal SNPs
			Arg6: n: int, sample size
			Arg7: list of strings, rsIDs
			Arg8: string, name of cluster (arbitrary)
			Arg9: numpy.array (2D), functional annotation matrix
			Arg10: numpy.array (1D), annotation effect size which will be estimated by EM
			Arg prior: string, independence_robust, "independence" or "gprior"
			Arg corr_thresh: int, excluding configurations with correlation > corr_thresh
			Arg v_scale = float, prior variance of the independence prior, recommended 0.05**2 (following FINEMAP, Benner et al 2016)
			Arg g: string, g-parameter of the g-prior, recommended g="BRIC" g=max(n,#SNPs**2), other options g="BIC" where g=n (Bayes Information Criterion), or  g="RIC" where g=#SNPs**2 (Risk Inflation Criterion) see Mixtures of g Priors for Bayesian Variable Selection Liang et al 2008
	'''
	nh_size = len(configs)
	configuration_size = numpy.array([len(configuration) for configuration in configs])

	# binomial prior without functional annotations
	if postgap.Globals.TYPE == 'binom':
		log_prior = calc_logbinom(configuration_size, kmax, len(z_scores))
	else:
		log_prior = calc_config_loglogis_priors(annotations, lambdas, configs, numpy.mean(
			numpy.exp(calc_logbinom(configuration_size, kmax, len(z_scores))), axis=0))

	log_BF = numpy.zeros(nh_size)
	i = 0

	for configuration in configs:
		if tuple(configuration) in score_cache:
			log_BF[i] = score_cache[tuple(configuration)]
			i = i+1
			continue

		z_tuple = numpy.take(z_scores, configuration)
		cor_tuple = numpy.take(cor_scores, configuration)
		cov_tuple = cov_matrix[numpy.ix_(configuration, configuration)]
		if numpy.max(numpy.absolute(numpy.tril(cov_tuple, k=-1))) > corr_thresh:
			score_cache[tuple(configuration)] = 0
			i = i+1
			continue
		tuple_size = len(z_tuple)
		if prior == "independence":
			v = numpy.eye(tuple_size) * v_scale
			log_BF[i] = calc_logBF(z_tuple, cov_tuple, v, n)
		elif prior == "independence_robust":
			v = numpy.eye(tuple_size) * v_scale
			log_BF[i] = calc_robust_logBF(
				z_tuple, cor_tuple, cov_tuple, v, n, eigen_thresh)
		elif prior == "gprior":
			log_BF[i] = calc_loggBF(z_tuple, cov_tuple, n, g=g)
		else:
			assert False, "%s is not one of independence or gprior" % (prior)
		score_cache[tuple(configuration)] = log_BF[i]
		i = i+1

	return OneDConfigurationSample(
		configurations=dict((tuple(config), index) for index, config in enumerate(configs)),
		posterior=(log_BF + log_prior),
		log_BF=log_BF,
		configuration_size=configuration_size,
		log_prior=log_prior,
		labels=labels,
		sample_label=sample_label
	)


def calc_logBF(z, cov, v, n):
	'''
			Compute Bayes Factors with assumption of independent variances
			Arg1: numpy.array (1D)
			Arg2: numpy.array (2D)
			Arg3: numpy.array diagonal matrix of prior variances
			Arg4: int sample size
			Returntype: numpy.array
	'''
	z = numpy.matrix(z)
	v = numpy.matrix(v)
	m = z.shape[1]
	try:
		coeff = 1. / math.sqrt(numpy.linalg.det((numpy.matrix(numpy.eye(m)) + n * numpy.matrix(v) * numpy.matrix(cov))))
	except ValueError:
		eigen_vec = numpy.linalg.eig(numpy.matrix(cov))[0]
		Q = numpy.linalg.eig(numpy.matrix(cov))[1]
		eigen_pos = eigen_vec
		eigen_pos[eigen_vec < 0] = 0
		cov = Q*numpy.diag(eigen_pos)*Q.transpose()
		coeff = 1. / math.sqrt(numpy.linalg.det((numpy.matrix(numpy.eye(m)) + n * numpy.matrix(v) * numpy.matrix(cov))))
	exponent = 0.5 * z * numpy.matrix(numpy.linalg.pinv(((n * v).I + cov), 0.0001)) * z.T
	return numpy.array(numpy.log(coeff) + exponent)[0][0]

def calc_robust_logBF(z, cor, cov, v, n, eigen_thresh):
	'''
			Compute Bayes Factors with assumption of independent variances, takes care of issues with negative definit correlation matrices
			%https://www.nag.co.uk/IndustryArticles/fixing-a-broken-correlation-matrix.pdf
			Arg1: numpy.array (1D)
			Arg2: numpy.array (2D)
			Arg3: numpy.array diagonal matrix of prior variances
			Arg4: int sample size
			Returntype: numpy.array
	'''
	z = numpy.matrix(z)
	cor = numpy.matrix(cor)
	v = numpy.matrix(v)
	m = z.shape[1]
	try:
		coeff = 1. / math.sqrt(numpy.linalg.det((numpy.matrix(numpy.eye(m)) + n * numpy.matrix(v) * numpy.matrix(cov))))
	except ValueError:
		eigen_vec = numpy.linalg.eig(numpy.matrix(cov))[0]
		Q = numpy.linalg.eig(numpy.matrix(cov))[1]
		eigen_pos = eigen_vec
		eigen_pos[eigen_vec < 0] = 0
		cov = Q*numpy.diag(eigen_pos)*Q.transpose()
		coeff = 1. / math.sqrt(numpy.linalg.det((numpy.matrix(numpy.eye(m)) + n * numpy.matrix(v) * numpy.matrix(cov))))
	exponent = 0.5 * z * \
		numpy.matrix(numpy.linalg.pinv(((n * v).I + cov), 0.0001)) * z.T
	log_BF = numpy.array(numpy.log(coeff) + exponent)[0][0]

	if(m > 1 and check_eigenvals(cor, cov, eigen_thresh) == False):
		return 0.
	else:
		return log_BF

def calc_loggBF(z, cov, n, g="BRIC"):
	'''
			Compute Bayes Factors with gprior method
			http://ftp.isds.duke.edu/WorkingPapers/05-12.pdf
			Arg1: numpy.array (1D)
			Arg2: numpy.array (2D)
			Arg3: int sample size
			Arg g: string, g-parameter of the g-prior, recommended g="BRIC" g=max(n,#SNPs**2), other options g="BIC" where g=n (Bayes Information Criterion), or  g="RIC" where g=#SNPs**2 (Risk Inflation Criterion) see Mixtures of g Priors for Bayesian Variable Selection Liang et al 2008
			Returntype: numpy.array
	'''
	if(g == "BRIC"):
		gp = numpy.max((len(z), n))
	if(g == "BIC"):
		gp = n
	if(g == "RIC"):
		gp = len(z)
	z = numpy.matrix(z)
	cov = numpy.matrix(cov)
	m = numpy.float(len(z))
	gp = numpy.float(gp)
	pinv = numpy.matrix(numpy.linalg.pinv(cov, 0.0001))
	coeff = (1 + gp)**(-m / 2)
	exponent = 0.5 * numpy.divide(gp, (gp + 1)) * z * pinv * z.T
	return numpy.array((math.log(coeff) + exponent))[0][0]

def calc_logbinom(subset_size, k, m):
	'''
			creates a binomial prior for a given subset_size under the assumption of k causal variants 
			excluding k=0 (that is no causal variant)
			Arg1: numpy.array (1D)
			Arg2: int, maximum number of causal variants
			Arg3: int, number of SNPs (len(z))
			Returntype: float
	'''
	p = float(1) / m

	if k == 1:
		return subset_size*numpy.log(p) + (m - subset_size)*numpy.log(1-p)
	else:
		logsum = 1*numpy.log(p) + (m - 1)*numpy.log(1-p)
		for i in range(2, k+1):
		    logsum = sumlog(logsum, (i*numpy.log(p)+(m-1)*numpy.log(1-p)))

		return (subset_size*numpy.log(p) + (m - subset_size)*numpy.log(1-p)) - logsum

def merge_samples(results, results_nh, labels, sample_label):
	'''
		Return merged OneDConfigurationSample
		Arg1: [ OneDConfigurationSample ]
		Arg2: [ OneDConfigurationSample ]
		Returntype: OneDConfigurationSample
	'''
	new_cfs = set(key for key in results_nh.configurations if key not in results.configurations)
	results.configurations.update(dict((cf, len(results.configurations) + i) for i,cf in enumerate(new_cfs)))
	to_add = [results_nh.configurations[cf] for cf in new_cfs]
	return OneDConfigurationSample(
		configurations = results.configurations,
		posterior = numpy.append(results.posterior, results_nh.posterior[to_add]),
		log_BF = numpy.append(results.log_BF, results_nh.log_BF[to_add]),
		configuration_size = numpy.append(results.configuration_size, results_nh.configuration_size[to_add]),
		log_prior = numpy.append(results.log_prior, results_nh.log_prior[to_add]),
		labels = labels,
		sample_label = sample_label
	)


def check_eigenvals(cor, Sigma, eigen_thresh=0.1):
	''' input 
	# z: vector of z-scores (aligned to the same reference as Sigma)
	# Sigma: correlation matrix
	'''

	assert cor.shape[1] == Sigma.shape[1], 'Covariance matrix has %i rows, %i expected' % (
		Sigma.shape[1], cor.shape[1])

	dim = cor.shape[1]
	S_tanh = numpy.matrix(numpy.zeros((dim+1, dim+1)))
	S_tanh[0, 0] = 1
	S_tanh[1:dim+1, 1:dim+1] = Sigma
	S_tanh[0, 1:dim+1] = cor
	S_tanh[1:dim+1, 0] = numpy.transpose(cor)

	return(numpy.min(numpy.linalg.eigvals(S_tanh)) > eigen_thresh)


def calc_approx_v(maf, sampN):
	'''
		compute approximate prior
		Arg1: float, Minor Allele Frequency (MAF) at SNP
		Arg2: int, sample size at SNP
		Returntype: float, approximate prior, v
	'''
	if maf != 0.:
		approx_v = 1.0 * maf * (1.0 - maf) * sampN
		approx_v = 1.0 / approx_v
	else:
		# approx_v = .0001
		approx_v = 10.0
	return approx_v


def calc_approx_v_cc(maf, sampN_case, sampN_control):
	'''
		compute approximate prior for case-control studies
		Arg1: float, Minor Allele Frequency (MAF) at SNP
		Arg2: int, sample size for case at SNP
		Arg3: int, sample size for control SNP
	'''
	if maf > 0.5:
		maf = 1 - maf
	if maf != 0.:
		num = sampN_case + sampN_control
		tmp1 = 2*maf*(1-maf) + 4*maf*maf
		interior = 2*maf*(1-maf) + 2*maf*maf
		tmp2 = interior * interior
		denom = sampN_case * sampN_control * (tmp1 - tmp2)
		return num / denom
	else:
		return .0001

def calc_logBF_ind(v, w, z_score):
	'''
		compute log Bayes Factor for one specified variance of prior(W)
		Arg1: float, approximate v prior at SNP
		Arg2: float, predefined variance of prior, w
		Arg3: float, z_score at SNP
		Returntype: float, log Bayes factor in the case of single prior, w
	'''
	r = w / (v + w)
	return -numpy.log(numpy.sqrt(1.0 - r)) - (z_score * z_score * r / 2.0)


def sumlog(logx, logy):
	'''
		compute log(x + y) = log( exp( log(x) ) + exp( log(y) ) )
	'''
	if logx > logy:
		return logx + numpy.log(1 + numpy.exp(logy - logx))
	else:
		return logy + numpy.log(1 + numpy.exp(logx - logy))


def calc_logBF_ML(v, W, z_score):
	'''
		compute averaged prior
		Arg1: float, approximate v prior at SNP
		Arg2: numpy.array (1D), prior variances
		Arg3: float, z_score at SNP
		Returntype: float, averaged snp log Bayes factor over W vector
	'''
	logBF_SNP = calc_logBF_ind(v, W[0], z_score)
	for i in range(1, len(W)):
	   logBF_SNP = sumlog(logBF_SNP, calc_logBF_ind(v, W[i], z_score))
	return logBF_SNP - numpy.log(len(W))

def calc_logBF(z, cov, v, n):
	'''
			Compute Bayes Factors with assumption of independent variances
			Arg1: numpy.array (1D)
			Arg2: numpy.array (2D)
			Arg3: numpy.array diagonal matrix of prior variances
			Arg4: int sample size
			Returntype: numpy.array
	'''
	z = numpy.matrix(z)
	v = numpy.matrix(v)
	m = z.shape[1]
	try:
		coeff = 1. / math.sqrt(numpy.linalg.det((numpy.matrix(numpy.eye(m)) + n * numpy.matrix(v) * numpy.matrix(cov))))
	except ValueError:
		eigen_vec = numpy.linalg.eig(numpy.matrix(cov))[0]
		Q = numpy.linalg.eig(numpy.matrix(cov))[1]
		eigen_pos = eigen_vec
		eigen_pos[eigen_vec < 0] = 0
		cov = Q*numpy.diag(eigen_pos)*Q.transpose()
		coeff = 1. / math.sqrt(numpy.linalg.det((numpy.matrix(numpy.eye(m)) + n * numpy.matrix(v) * numpy.matrix(cov))))
	exponent = 0.5 * z * numpy.matrix(numpy.linalg.pinv(((n * v).I + cov), 0.0001)) * z.T
	return numpy.array(numpy.log(coeff) + exponent)[0][0]


def set_prior(pi, annot, lambdas):
	'''
		compute logistic prior at SNP
		we assume intecept is based on probability of pi which can be predefined by user
		So, we don't estimate this parameter in the model but we put this value in the model
		Arg1: float, a predefined probability for intercept where there is no annotation at SNP 
		Arg2: numpy.array (1D), annotation information at SNP
		Arg3: numpy.array (1D), annotation effect sizes
		Returntype: float, logistic prior probability at SNP
	'''
	assert pi > 0 and pi < 1
	annot = annot.transpose()
	logitprior = numpy.log(pi) - numpy.log(1 - pi) + sum(annot * lambdas)
	return 1.0 / (1.0 + numpy.exp(-logitprior))

def set_prior_prob(pi, annot, lambdas):
	'''
		compute logistic prior at SNP
		we assume intecept is based on probability of pi which can be predefined by user
		So, we don't estimate this parameter in the model but we put this value in the model
		Arg1: float, a predefined probability for intercept where there is no annotation at SNP 
		Arg2: numpy.array (1D), annotation information at SNP
		Arg3: numpy.array (1D), annotation effect sizes
		Returntype: float, logistic prior probability at SNP
	'''
	annot = annot.transpose()
	logitprior = numpy.log(pi) - numpy.log(1 - pi) + sum(annot * lambdas)
	return 1.0 / (1.0 + numpy.exp(-logitprior))

def compute_eqtl_lambdas(cluster, z_scores, kmax, W=[0.01, 0.1, 0.5], pi=0.01):
	'''
		The function for leanring F.A parameters in eQTL
		Arg1 GWAS_Cluster
		Arg2: Numpy Array
		Arg3: [float]
		Arg4: float
		Returntype: numpy.array
	'''
	initial_lambdas = [0.0] * cluster.annotations.shape[0]
	out = calc_logbinom(
		cluster.ld_matrix.shape[0], len(cluster.z_scores))
	pi = numpy.exp(out)
	MAFs = map(float, cluster.mafs)
	if postgap.Globals.TYPE == 'binom' or postgap.Globals.TYPE == 'EM':
		return initial_lambdas
	else:
		sample_size = cluster.gwas_snps[0].evidence[0].sample_size
		approx_v = [calc_approx_v(maf, sample_size) for maf in MAFs]
		logBFs = [calc_logBF_ML(approx_v[i], W, z_scores[i])
				 for i in range(len(z_scores))]
		mat_annot = cluster.annotations.T

		def f_llk(lambdas_):
			priors = [set_prior(pi, annot, lambdas_) for annot in mat_annot]
			lsum = 0
			for i in range(len(priors)):
				if priors[i] < 1.2683388472255632e-301:
					priors[i] = 1.2683388472255632e-301
				if priors[i] > 0.9999999999999999:
					priors[i] = 0.9999999999999999
				lsum = lsum + \
					sumlog(logBFs[i] + numpy.log(priors[i]),
						   numpy.log(1 - priors[i]))
			return 1/lsum
		bounds_list = [(-30, 30) for i in range(cluster.annotations.shape[0])]
		result = scipy.optimize.minimize(f_llk, initial_lambdas, method='L-BFGS-B', bounds=bounds_list)
		if result.success:
			return 1 / (1 + numpy.exp(-result.x))
		else:
			return initial_lambdas

def compute_gwas_lambdas(cluster, kmax, W=[0.01, 0.1, 0.5]):
	'''
		compute MLE(Maximum Likelihood Estimate) of annotation effect size
		build clusters with this MLE
		Arg1: GWAS_Cluster
		Arg2: variance of prior which will be averaged over [0.01, 0.1, 0.5]
		Returntype: GWAS_Cluster
		'''
	pi = numpy.exp(calc_logbinom(
		cluster.ld_matrix.shape[0], kmax, len(cluster.z_scores)))
	initial_lambdas = [0.0] * cluster.annotations.shape[0]
	MAFs = map(float, cluster.mafs)
	if postgap.Globals.TYPE == 'binom' or postgap.Globals.TYPE == 'EM':
		lambdas = initial_lambdas
	else:
		sample_size = cluster.gwas_snps[0].evidence[0].sample_size
		# To DO: Setting option for quantative and qualitative trait (Default for GWAS is qualitative and for eQTL is quantitative)
		approx_v = [calc_approx_v(maf, sample_size) for maf in MAFs]
		# To DO: sampN_case / sampN_control setting

		logBFs = [calc_logBF_ML(approx_v[i], W, cluster.z_scores[i])
				  for i in range(len(cluster.z_scores))]
		mat_annot = cluster.annotations.T
		assert len(mat_annot) == len(logBFs), "Matrices of discrepant sizes: annotations.T = %i; logBFs = %i, cluster.annotations.shape = %s, ld.matrix = %s" % (len(cluster.annotations.T), len(logBFs), cluster.annotations.shape, cluster.ld_matrix.shape)
		min_val = 1.2683388472255632e-301
		max_val = 0.9999999999999999

		def f_llk(lambdas_):
			priors = [set_prior(pi, annot, lambdas_) for annot in mat_annot]
			lsum = 0
			for i in range(len(priors)):
				if priors[i] < min_val:
					priors[i] = min_val
				if priors[i] > max_val:
					priors[i] = max_val
				lsum = lsum + \
					sumlog(logBFs[i] + numpy.log(priors[i]),
						   numpy.log(1 - priors[i]))
			return 1/lsum

		bounds_list = [(-30, 30) for i in range(cluster.annotations.shape[0])]
		result = scipy.optimize.minimize(
			f_llk, initial_lambdas, method='L-BFGS-B', bounds=bounds_list)
		if result.success:
			lambdas = 1.0 / (1.0 + numpy.exp(-result.x))
		else:
			lambdas = initial_lambdas

	return GWAS_Cluster(cluster.gwas_snps,
			 cluster.ld_snps,
			 cluster.ld_matrix,
			 cluster.z_scores,
			 cluster.betas,
			 MAFs,
			 cluster.annotations,
			 None,
			 lambdas)


def calc_config_loglogis_priors(annotations, lambdas, configurations, pi_EM=0.001):
	'''
		compute log logistic prior at config
		Arg1: numpy.array (2D), functional annotations
		Arg2: numpy.array (1D), annotation effect sizes
		Arg3: array of arrays
		Arg pi_EM: predefined parameter for setting intercept in logistic prior during EM calculation
		Returntype: Numpy.array
	'''
	return numpy.array([calc_config_loglogis_prior(annotations, lambdas, configuration, pi_EM) for configuration in configurations])

def calc_config_loglogis_prior(annotations, lambdas, configuration, pi_EM=0.001):
	'''
		compute log logistic prior at config
		Arg1: numpy.array (2D), functional annotations
		Arg2: numpy.array (1D), annotation effect sizes
		Arg3: array 
		Arg pi_EM: predefined parameter for setting intercept in logistic prior during EM calculation
		Returntype: Numpy.array
	'''
	logprior = 0
	for snp in range(annotations.shape[1]):
		set_prior_out = set_prior_prob(pi_EM, annotations[:, snp], lambdas)
		if snp in configuration:
			if set_prior_out < 1e-301:
				logprior = logprior + numpy.log(1e-301)
			else:
				logprior = logprior + numpy.log(set_prior_out)
		else:
			if set_prior_out > (1-1e-16):
				logprior = logprior + numpy.log(1 - (1-1e-16))
			else:
				logprior = logprior + numpy.log(1 - set_prior_out)
	
	if logprior < -700:
		return -700
	else:
		return logprior
