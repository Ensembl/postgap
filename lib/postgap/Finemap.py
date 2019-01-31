# functions for one-dimensional and two-dimensional fine-mapping
# shotgun stochastic search
# 11 april 2017


import numpy
import math
import itertools as it
import operator
import random
import collections

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
		'''
			Return OneDConfigurationSample where the configuration posteriors are reduced to marginals
			on the positions
			Arg2: OneDConfigurationSample
			Returntype: OneDConfigurationSample
		'''
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
		'''
			INPUT----
			OneDConfigurationSample self
			Integer kmax
			Integer singleton_count
			OUTPUT------
			Dictionary key range[1, kmax] value OneDConfigurationSample
		'''
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
			      	configuration =  input_keys[count]      
				#print(configuration)
				if  (len(configuration)<nc): 
	  				count += 1
				else: 
	  				index_tupel =  [config for config in it.combinations(configuration,nc)]
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
		'''
			Main function for fine-mapping of two traits
			Arg1: OneDConfigurationSample
			Arg2: OneDConfigurationSample
			Returntype: (float (the overall evidence for shared causal variation at a locus), TwoDConfigurationSample)
		'''
		
		sample1_labels = self.labels
		sample2_labels = sample2.labels
		
		translate_from_sample2_to_sample1 = self.create_translation_index(
			from_labels = sample2_labels,
			to_labels   = sample1_labels
		)
		
		translated_sample2_configurations = self.translate_configurations(
			translation_index = translate_from_sample2_to_sample1,
			configurations    = sample2.configurations
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
		indentation = "    "
		total_configurations = len(self.configurations)
		if total_configurations >  max_show_at_top + max_show_at_end:
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
			summary_lines.append(indentation +  "...")

		for configuration in sorted_configurations[-num_show_at_end:]:
			summary_lines.append(indentation + self.configuration_string(configuration))

		return "\n".join(summary_lines)

	def configuration_string(self, configuration):
		index_of_configuration = configurations[configuration]
		posterior = self.posterior[index_of_configuration]
		prior     = math.exp(self.log_prior[index_of_configuration])
		BF	= math.exp(log_BF[index_of_configuration])
		snps      = ', '.join([self.labels[position] for position in configuration])
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

	def string(self, max_show_at_top = 5, max_show_at_end = 2):
		indentation = "    "
		total_configurations = len(self.configurations)
		if total_configurations >  max_show_at_top + max_show_at_end:
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
			summary_lines.append(indentation +  "...")

		for configuration in sorted_configurations[-num_show_at_end:]:
			summary_lines.append(indentation + self.configuration_string(configuration))

		summary_lines.append("The coloc evidence is: {:2e}".format(self.coloc_evidence))
		
		return '\n'.join(summary_lines)

	def configuration_string(self, configuration):
		snp_labels = ", ".join([self.labels[position] for position in configuration])
		index = self.configurations[configuration]
		posterior  = self.posterior[index]
		posterior1 = self.posterior1[index]
		posterior2 = self.posterior2[index]
		return "The snp configuration ({}) had posterior probabilities of {:.2e} in {} and {:.2e} in {}. The joint posterior probability is {:.2e}.".format(
			snp_labels,
			float(posterior1),
			self.sample_1_label,
			float(posterior2),
			self.sample_2_label,
			float(posterior))

def finemap(z_scores, beta_scores, cov_matrix, n, labels, sample_label, mafs = None, annotations=None, kstart=1, kmax=5, corr_thresh=0.9, max_iter=100000, output="configuration", prior="independence_robust", v_scale=0.0025, g="BRIC", eigen_thresh=0.1, verbose=False):
	'''
		Main function for fine-mapping using stochastic search for one trait #
		Arg1: z_scores: numpy.array
		Arg2: beta_scores: numpy.array
		Arg3: cov_matrix numpy.array, correlation-structure (TwoD) array
		Arg4: n: int, sample size
		Arg5: list of strings, snp rsIDs
		Arg6: name of cluster (arbitrary) for easier reporting
		Arg kstart: int, full exploration of sets with #kstart causal variants
		Arg kmax: int, maximum number of causal variants
		Arg corr_thresh: int, excluding configurations with correlation > corr_thresh
		Arg max_iter: int, iterations of stochastic search
                Arg output: int, "configuration" or "marginal"
		Arg prior= string ("independence", "independence_robust" or "gprior") independence_robust does filter out configurations with possible mis-matches between the correlation matrix and effect size direction
		Arg v_scale = float, prior variance of the independence prior, recommended 0.05 **2 (following FINEMAP, Benner et al 2016)
		Arg g: string, g-parameter of the g-prior, recommended g="BRIC" g=max(n,#SNPs**2), other options g="BIC" where g=n (Bayes Information Criterion), or  g="RIC" where g=#SNPs**2 (Risk Inflation Criterion) see Mixtures of g Priors for Bayesian Variable Selection Liang et al 2008
		Arg verbose = False: boolean, print the progress of the stochastic search
		Returntype: OneDConfigurationSample
	'''

	# Test inputs
	assert len(z_scores) == cov_matrix.shape[0], 'Covariance matrix has %i rows, %i expcted' % (cov_matrix.shape[0], len(z_scores))
	assert len(z_scores) == cov_matrix.shape[1], 'Covariance matrix has %i columns, %i expcted' % (cov_matrix.shape[0], len(z_scores))
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
		configs  = configurations, 
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
			results_nh = compare_neighborhood(new_configs, z_scores, cor_scores, cov_matrix, kmax, n, score_cache, prior, corr_thresh,  v_scale=v_scale, g=g)

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

		res_out = merge_samples(result_list).normalise_posteriors()

		if output == "configuration":
			return res_out
		elif output == "marginal":
			return res_out.marginals(len(z_scores))
		else:
			assert False, "%s unkown" % {output}

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

def compare_neighborhood(configs, z_scores,  cor_scores, cov_matrix, kmax, n, score_cache, labels, sample_label, prior="independence_robust", corr_thresh=0.9, v_scale=0.0025, g="BRIC", eigen_thresh=0.1):
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
		Arg prior: string, independence_robust, "independence" or "gprior"
		Arg corr_thresh: int, excluding configurations with correlation > corr_thresh
		Arg v_scale = float, prior variance of the independence prior, recommended 0.05**2 (following FINEMAP, Benner et al 2016)
		Arg g: string, g-parameter of the g-prior, recommended g="BRIC" g=max(n,#SNPs**2), other options g="BIC" where g=n (Bayes Information Criterion), or  g="RIC" where g=#SNPs**2 (Risk Inflation Criterion) see Mixtures of g Priors for Bayesian Variable Selection Liang et al 2008
	'''
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
	'''
		Compute Bayes Factors with gprior method
		http://ftp.isds.duke.edu/WorkingPapers/05-12.pdf
		Arg1: numpy.array (1D)
		Arg2: numpy.array (2D)
		Arg3: int sample size
		Arg g: string, g-parameter of the g-prior, recommended g="BRIC" g=max(n,#SNPs**2), other options g="BIC" where g=n (Bayes Information Criterion), or  g="RIC" where g=#SNPs**2 (Risk Inflation Criterion) see Mixtures of g Priors for Bayesian Variable Selection Liang et al 2008
		Returntype: numpy.array
	'''
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
	'''
		creates a binomial prior for a given subset_size under the assumption of k causal variants 
		excluding k=0 (that is no causal variant)
		Arg1: numpy.array (1D)
		Arg2: int, maximum number of causal variants
		Arg3: int, number of SNPs (len(z))
		Returntype: float
	'''
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

def merge_samples(samples):
	'''
		Return merged OneDConfigurationSample
		Arg1: [ OneDConfigurationSample ]
		Returntype: OneDConfigurationSample
	'''
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
			log_prior = log_prior
		)

def check_eigenvals(cor,Sigma,eigen_thresh=0.1):
	''' input 
	# z: vector of z-scores (aligned to the same reference as Sigma)
	# Sigma: correlation matrix
	'''

	assert cor.shape[1] == Sigma.shape[1], 'Covariance matrix has %i rows, %i expected' % (Sigma.shape[1], cor.shape[1])

	dim=cor.shape[1]
	S_tanh = numpy.matrix(numpy.zeros((dim+1, dim+1)))
	S_tanh[0,0] = 1
	S_tanh[1:dim+1,1:dim+1] = Sigma
	S_tanh[0,1:dim+1] = cor
	S_tanh[1:dim+1,0] = numpy.transpose(cor)

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
        approx_v = .001
    return approx_v

def calc_logBF_ind(v, w, z_score):
    '''
        compute log Bayes Factor for one specified variance of prior(W)
        Arg1: float, approximate v prior at SNP
        Arg2: float, predefined variance of prior, w
        Arg3: float, z_score at SNP
        Returntype: float, log Bayes factor in the case of single prior, w
    '''
    r = w / (v + w)
    logBF_ind = numpy.log(numpy.sqrt(1.0 - r)) + (z_score * z_score * r / 2.0)
    return logBF_ind

def sumlog(logx, logy):
    '''
        compute log(x + y) = log( exp( log(x) ) + exp( log(y) ) )
    '''
    if logx > logy:
        return logx + numpy.log(1 + numpy.exp(logy - logx))
    else:
        return logy + numpy.log(1 + numpy.exp(logx - logy))

def calc_logBF(v, W, z_score):
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
    logBF_SNP = logBF_SNP - numpy.log(len(W))
    return logBF_SNP

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
    logitprior = numpy.log(pi) - numpy.log(1 - pi)
    logitprior = logitprior + sum( annot * lambdas )
#     for i in range(len(annot)):
#         if annot[i]:
#             logitprior = logitprior + lambdas[i]
    prior = 1.0 / (1.0 + numpy.exp(-logitprior))
    return prior

GWAS_Cluster_prototype_with_lambdas = collections.namedtuple(
        'GWAS_Cluster',
        [
                'gwas_snps',
                'ld_snps',
                'ld_matrix',
                'z_scores',
                'betas',
                'mafs',
                'annotations',
                'gwas_configuration_posteriors',
                'lambdas'
        ]
)
class GWAS_Cluster_with_lambdas(GWAS_Cluster_prototype_with_lambdas):
    pass


def mk_modified_clusters(p_cluster, W = [0.01, 0.1, 0.5], pi = 0.001):
    '''
        compute MLE(Maximum Likelihood Estimate) of annotation effect size
        build clusters with this MLE
        Arg1: GWAS_Cluster
        Returntype: GWAS_Cluster_with_lambdas
        Arg W: variance of prior which will be averaged over [0.01, 0.1, 0.5]
        Arg pi: predefined parameter for setting intercept in logistic prior
        '''
    MAFs       = map(float, p_cluster.mafs)
    sample_size= p_cluster.gwas_snps[0].evidence[0].sample_size
    approx_v   = [calc_approx_v(maf, sample_size)  for maf in MAFs]
    z_scores   = p_cluster.z_scores
    logBFs     = [calc_logBF(approx_v[i], W, z_scores[i]) for i in range(len(z_scores))]
    mat_annot  = p_cluster.annotations.T

    def f_llk(lambdas):
        priors = [set_prior(pi, annot, lambdas) for annot in mat_annot]

        lsum = 0
        for i in range(len(priors)):
            lsum = lsum + sumlog(logBFs[i] + numpy.log(priors[i]), numpy.log(1 - priors[i]))
        return -lsum

    initial_lambdas = [0.000001] * mat_annot.shape[1]
    result = optimize.minimize(f_llk, initial_lambdas, method='L-BFGS-B')
    if result.success:
        lambdas= 1 / (1 + numpy.exp(-result.x))
    else:
        lambdas= initial_lambdas

    return GWAS_Cluster_with_lambdas(   p_cluster.gwas_snps,
                                        p_cluster.ld_snps,
                                        p_cluster.ld_matrix,
                                        z_scores,
                                        p_cluster.betas,
                                        MAFs,
                                        p_cluster.annotations,
                                        None,
                                        lambdas)


