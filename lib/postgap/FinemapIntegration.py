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
# Hack to lazy load scipy only if required
from postgap.DataModel import *
import cPickle as pickle
import collections
import postgap.Globals
import numpy
import os.path
scipy = None


def compute_gwas_posteriors(cluster_associations, populations):
	"""
			Compute posterior of GWAS causality across all clusters
			Arg1: [(GWAS_Cluster, GeneSNP_Association)]
			Arg2: dict(string => float)
			Arg3: binom/ ML/ EM/ ML_EM
			Returntype: [(GWAS_Cluster, GeneSNP_Associations)]
	"""
	# prepped_clusters = [(prepare_cluster_for_finemap(cluster, associations, populations), associations) for cluster, associations in cluster_associations]
	prepped_clusters = []
	for cluster, associations in cluster_associations:
		if (len(cluster.ld_snps) < 3):
			continue
		prepped_clusters.append((prepare_cluster_for_finemap(
			cluster, associations, populations), associations))
		# prepped_cluster = prepare_cluster_for_finemap(
		#	 cluster, associations, populations)
		# prepped_cluster_rsIDs = [
		#	 ld_snp.rsID for ld_snp in prepped_cluster.ld_snps]
		# filtered_associations = [
		#	 association for association in associations if association.snp.rsID in prepped_cluster_rsIDs]
		# prepped_clusters.append((prepped_cluster, filtered_associations))

	# MLE calculation
	modified_clusters = [(postgap.Finemap.mk_modified_clusters(
		cluster), associations) for cluster, associations in prepped_clusters]
	pickle.dump(modified_clusters, open('prepped_clusters_' +
										os.path.basename(postgap.Globals.GWAS_SUMMARY_STATS_FILE)+'.pkl', "w"))
	return [(finemap_gwas_cluster(cluster), associations) for cluster, associations in modified_clusters]


def prepare_cluster_for_finemap(cluster, associations, population, tissue_weights=['Whole_Blood']):
	'''

			Enriches GWAS clusters with z-scores, betas, MAFs and annotations
			Arg1: GWAS_Cluster
			Arg2: [GeneSNP_Association]
			Arg2: dict(string => float)
			Returntype: GWAS_Cluster

	'''

	if len(cluster.ld_snps) == len(cluster.gwas_snps):
		ld_snps, ld_matrix, z_scores, betas = compute_ld_matrix(
			cluster, population)

	elif postgap.Globals.GWAS_SUMMARY_STATS_FILE is not None:
		ld_snps, ld_matrix, z_scores, betas = extract_z_scores_from_file(
			cluster, population)

	else:
		ld_snps, ld_matrix, z_scores, betas = impute_z_scores(cluster, population)
	
	# if ld_snps has less than 10 element, it is less informative
	# TOTO: TODO assert len(ld_snps) >= 10, sample_label + ' has less than 10 ld_snps in the cluster'
	
	mafs = extract_snp_mafs(ld_snps, associations, population.lower())
	annotations = (extract_snp_annotations(ld_snps, associations) > 0.).astype('float')

	assert len(ld_snps) == ld_matrix.shape[0]
	assert len(ld_snps) == ld_matrix.shape[1]
	return GWAS_Cluster(cluster.gwas_snps, ld_snps, ld_matrix, z_scores, betas, mafs, annotations, None)
	
def extract_snp_mafs(ld_snps, associations, populations):
	"""
			Produce vector of mafs for the SNPs provided
			Arg1: [SNP]
			Arg2: [GeneSNP_Association]
			Arg3: String
			Returntype: Numpy vector
	"""
	maf_hash = collections.defaultdict(int)
	for association in associations:
		for evidence in association.regulatory_evidence:
			if evidence.info is not None:
				if evidence.info['MAFs'] is not None:
					if evidence.source == 'VEP_reg' and populations in evidence.info['MAFs']:
						maf_hash[association.snp.rsID] = evidence.info['MAFs'][populations]

	return numpy.array([float(maf_hash[snp.rsID]) for snp in ld_snps])


def extract_snp_annotations(ld_snps, associations):
	"""
			Produce array of annotation for the SNPs provided 
			Arg1: [SNP]
			Arg2: [GeneSNP_Association]
			Returntype: Numpy 2D array
	"""
	annotation_hash = collections.defaultdict(lambda: collections.defaultdict(float))

	for association in associations:
		for evidence in association.cisregulatory_evidence + association.regulatory_evidence:
			if evidence.source in ['GTEx']:
				continue
			else:
				annotation_hash[evidence.source][evidence.snp.rsID] = evidence.score

	return numpy.array([[annotation_hash[annotation][snp.rsID] for snp in ld_snps] for annotation in sorted(annotation_hash.keys())])

def compute_ld_matrix(cluster, population):
	'''
			Computes LD matrix, re-orders SNPs zccordingly, and extract Z-scores from 
			Cluster data.
			Arg1: Cluster
			Arg2: population (string)
			Returntype: [SNP], numpy.matrix (square LD matrix), numpy.matrix (Z-score vector)
	'''
	# Compute ld_matrix
	ld_snp_ids, ld_matrix = postgap.LD.get_pairwise_ld(cluster.ld_snps, population)

	# Update list of LD SNPs
	ld_snp_hash = dict((ld_snp.rsID, ld_snp) for index, ld_snp in enumerate(cluster.ld_snps))
	ld_snps = [ld_snp_hash[rsID] for rsID in ld_snp_ids]

	# Aggregate z_scores into a single vector
	z_scores = [0] * len(ld_snps)
	betas = [0] * len(ld_snps)

	gwas_snp_hash = dict((gwas_snp.snp.rsID, gwas_snp) for gwas_snp in cluster.gwas_snps)
	for index, ld_snp in enumerate(ld_snps):
		z_scores[index] = gwas_snp_hash[ld_snp.rsID].z_score
		betas[index] = gwas_snp_hash[ld_snp.rsID].beta

	assert len(ld_snps) == ld_matrix.shape[0]
	assert len(ld_snps) == ld_matrix.shape[1]
	return ld_snps, ld_matrix, z_scores, betas


def extract_z_scores_from_file(cluster, population):
	'''
			Extracts Z-scores from summary stats file, computes LD matrix
			TODO: It is inefficient torun through a 1GB file once per locus, this should be done only once
			Arg1: Cluster
			Returntype: [SNP], numpy.matrix (square LD matrix), numpy.matrix (Z-score vector)
	'''
	ld_snp_hash = dict((ld_snp.rsID, ld_snp)
					   for index, ld_snp in enumerate(cluster.ld_snps))

	## Search ld snps in the GWAS summary stats file, and drop ld snps that cannot be found in the GWAS summary stats file from ld_snps
	proper_gwas_cluster = postgap.GWAS.GWAS_File().create_gwas_cluster_with_pvalues_from_file(gwas_cluster=cluster, gwas_data_file=postgap.Globals.GWAS_SUMMARY_STATS_FILE)

	## Extract z_scores for all the ld snps that can be found in the GWAS summary stats file
	ld_snp_results = dict((ld_snp.snp, (ld_snp.pvalue, ld_snp.beta)) for index, ld_snp in enumerate(proper_gwas_cluster.ld_snps))
	
	# Select all the found ld snps from ld_snps
	found_ld_snps = [ld_snp_hash[rsID] for rsID in ld_snp_results]

	# Compute ld_matrix
	ld_snp_ids, ld_matrix = postgap.LD.get_pairwise_ld(
		found_ld_snps, population)

	# Update list of LD SNPs
	ld_snps = [ld_snp_hash[rsID] for rsID in ld_snp_ids]
	z_scores = [z_score_from_pvalue(ld_snp_results[rsID][0], ld_snp_results[rsID][1]) for rsID in ld_snp_ids]
	betas = [ld_snp_results[rsID][1] for rsID in ld_snp_ids]

	assert len(ld_snps) == ld_matrix.shape[0]
	assert len(ld_snps) == ld_matrix.shape[1]
	return ld_snps, ld_matrix, z_scores, betas


def impute_z_scores(cluster, population):
	'''
			Imputes Z-scores from available data, computes LD matrix
			Arg1: Cluster
			Returntype: [SNP], numpy.matrix (square LD matrix), numpy.matrix (Z-score vector)
	'''
	# Compute ld_matrix
	ld_snp_ids, ld_matrix = postgap.LD.get_pairwise_ld(cluster.ld_snps, population)

	# Update list of LD SNPs
	ld_snp_hash = dict((ld_snp.rsID, ld_snp) for index, ld_snp in enumerate(cluster.ld_snps))
	ld_snps = [ld_snp_hash[rsID] for rsID in ld_snp_ids]

	# Determine which SNPs are missing values
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
	shrink_lambda = 0.1  # shrinkage factor, magic number
	ld_matrix_known_shrink = shrink_lambda * numpy.diag(numpy.ones(ld_matrix_known.shape[0])) + (1-shrink_lambda) * ld_matrix_known
	ld_matrix_k2m_shrink = (1-shrink_lambda) * ld_matrix_k2m
	z_shrink_imputed = numpy.dot(numpy.dot(ld_matrix_k2m_shrink, numpy.linalg.pinv(ld_matrix_known_shrink, 0.0001)), known_z_scores)
	beta_shrink_imputed = numpy.dot(numpy.dot(ld_matrix_k2m_shrink, numpy.linalg.pinv(ld_matrix_known_shrink, 0.0001)), known_betas)

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

	assert len(ld_snps) == ld_matrix.shape[0]
	assert len(ld_snps) == ld_matrix.shape[1]
	return ld_snps, ld_matrix, z_scores, betas

def finemap_gwas_cluster(cluster):
	'''

			Enriches GWAS clusters with z-scores and GWAS posteriors
			Arg1: GWAS_Cluster
			Returntype: GWAS_Cluster

	'''
	logging.info("Finemap GWAS Cluster")
	# Define experiment label (serves for debugging logs)
	chrom = cluster.ld_snps[0].chrom
	start = min(ld_snp.pos for ld_snp in cluster.ld_snps)
	end = max(ld_snp.pos for ld_snp in cluster.ld_snps)
	sample_label = 'GWAS_Cluster_%s:%i-%i' % (chrom, start, end)

	# Define LD SNP labels (serves for debugging logs)
	ld_snp_ids = [ld_snp.rsID for ld_snp in cluster.ld_snps]

	# Define sample size: mean of max for each SNP
	sample_sizes = map(lambda gwas_snp: max(
		gwas_association.sample_size for gwas_association in gwas_snp.evidence), cluster.gwas_snps)
	sample_size = sum(sample_sizes) / len(sample_sizes)

	# Compute posterior
	if postgap.Globals.TYPE == 'binom' or postgap.Globals.TYPE == 'ML':
		configuration_posteriors = postgap.Finemap.finemap_v1(
			z_scores=numpy.array(cluster.z_scores),
			beta_scores=numpy.array(cluster.betas),
			cov_matrix=cluster.ld_matrix,
			n=sample_size,
			labels=ld_snp_ids,
			sample_label=sample_label,
			lambdas=cluster.lambdas,
			mafs=cluster.mafs,
			annotations=cluster.annotations,
			kmax=postgap.Globals.kmax_gwas,
			isGWAS=True
		)
	elif postgap.Globals.TYPE == 'EM' or postgap.Globals.TYPE == 'ML_EM':
		configuration_posteriors = postgap.Finemap.finemap_v2(
			z_scores=numpy.array(cluster.z_scores),
			beta_scores=numpy.array(cluster.betas),
			cov_matrix=cluster.ld_matrix,
			n=sample_size,
			labels=ld_snp_ids,
			sample_label=sample_label,
			lambdas=cluster.lambdas,
			mafs=cluster.mafs,
			annotations=cluster.annotations,
			kmax=postgap.Globals.kmax_gwas,
			isGWAS=True
		)
	return GWAS_Cluster(cluster.gwas_snps, cluster.ld_snps, cluster.ld_matrix, cluster.z_scores, cluster.betas, cluster.mafs, cluster.annotations, configuration_posteriors, cluster.lambdas)


def compute_joint_posterior(cluster, gene_tissue_snp_eQTL_hash):
	"""
			Compute collocation posterior of gene expression and GWAS phenotype at the specified cluster and tissue
			Arg1: GWAS_Cluser
			Arg4: [GeneSNP_Association]
			Returntype: Hash of hashes: Gene => Tissue => (rsID|CLUSTER) => Float
	"""
	return dict((gene, compute_gene_joint_posterior(cluster, gene, gene_tissue_snp_eQTL_hash[gene], cluster.gwas_configuration_posteriors, cluster.mafs, cluster.annotations)) for gene in gene_tissue_snp_eQTL_hash)


def compute_gene_joint_posterior(cluster, gene, tissue_snp_eQTL_hash, gwas_configuration_posteriors, mafs, annotations):
	"""
			Compute collocation posterior of gene expression and GWAS phenotype at the specified cluster and tissue
			Arg1: GWAS_Cluster
			Arg2: Gene
			Arg3: Hash of hashes: Tissue => SNP => (Float, Float)
			Arg4: Hash of hashes: configuration => posterior
			Arg5: Numpy Vector
			Arg6: Numpy 2D Array
			Returntype: Hash of hashes: Tissue => Float
	"""
	assert len(cluster.ld_snps) == cluster.ld_matrix.shape[0], (len(cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
	assert len(cluster.ld_snps) == cluster.ld_matrix.shape[1], (len(cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])

	return dict((tissue, compute_gene_tissue_joint_posterior(cluster, gene, tissue, tissue_snp_eQTL_hash[tissue], gwas_configuration_posteriors, mafs, annotations)) for tissue in tissue_snp_eQTL_hash)


def compute_gene_tissue_joint_posterior(cluster, gene, tissue, eQTL_snp_hash, gwas_configuration_posteriors, mafs, annotations):
	"""
			Compute posterior of gene expression regulation at the specified cluster and tissue
			Arg1: GWAS_Cluster
			Arg2: Tissue (string)
			Arg3: Gene
			Arg4: Hash string (rsID) => (Float (z-score), Float (beta))
			Arg5: Hash of hashes: configuration => posterior
			Returntype: Numpy Array (GWAS PIP values), Numpy Array (Coloc evidence)
	"""
	# eQTL posteriors
	eQTL_configuration_posteriors = compute_eqtl_posteriors(cluster, tissue, gene, eQTL_snp_hash, mafs, annotations)

	# Joint posterior
	joint_out = eQTL_configuration_posteriors.joint_posterior(gwas_configuration_posteriors)

	return joint_out[2], joint_out[3]


def compute_eqtl_posteriors(cluster, tissue, gene, eQTL_snp_hash, mafs, annotations):
	"""
			Compute posterior of gene expression regulation at the specified cluster and tissue
			Arg1: GWAS_Cluster
			Arg2: Tissue (string)
			Arg3: Gene
			Arg4: Hash string (rsID) => (Float (z-score), Float (beta))
			Arg5: Numpy Vector
			Arg6: Numpy 2D Array
			Returntype: Float
	"""
	assert len(cluster.ld_snps) == cluster.ld_matrix.shape[0], (len(cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
	assert len(cluster.ld_snps) == cluster.ld_matrix.shape[1], (len(cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
	# Determine which SNPs are missing values
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
		shrink_lambda = 0.1  # shrinkage factor, magic number
		ld_matrix_known_shrink = shrink_lambda * numpy.diag(numpy.ones(ld_matrix_known.shape[0])) + (1-shrink_lambda) * ld_matrix_known
		assert ld_matrix_known_shrink.size > 0, (missing_indices, ld_matrix, ld_matrix_known)
		ld_matrix_k2m_shrink = (1-shrink_lambda) * ld_matrix_k2m
		z_shrink_imputed = numpy.dot(numpy.dot(ld_matrix_k2m_shrink, numpy.linalg.pinv(ld_matrix_known_shrink, 0.0001)), known_z_scores)
		beta_shrink_imputed = numpy.dot(numpy.dot(ld_matrix_k2m_shrink, numpy.linalg.pinv(ld_matrix_known_shrink, 0.0001)), known_betas)

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

	# Define experiment label (serves for debugging logs)
	chrom = cluster.ld_snps[0].chrom
	start = min(ld_snp.pos for ld_snp in cluster.ld_snps)
	end = max(ld_snp.pos for ld_snp in cluster.ld_snps)
	sample_label = 'eQTL_Cluster_%s:%i-%i_%s' % (chrom, start, end, gene)

	# Learn F.A parameters in eQTL 
	lambdas = postgap.Finemap.compute_eqtl_lambdas(cluster, numpy.array(z_scores))

	# Compute posterior
	logging.debug("Finemap eQTL Cluster")
	return postgap.Finemap.finemap_v1(
		z_scores=numpy.array(z_scores),
		beta_scores=numpy.array(betas),
		cov_matrix=cluster.ld_matrix,
		n=500,  # TODO extract eQTL sample sizes
		labels=[ld_snp.rsID for ld_snp in cluster.ld_snps],
		sample_label=sample_label,
		lambdas=lambdas,
		mafs=mafs,
		annotations=annotations,
 		kstart       = postgap.Globals.KSTART,
                kmax         = postgap.Globals.KMAX,
		isGWAS=False
	)

	## Joint posterior
	sum_posteriors, config_sample = eQTL_configuration_posteriors.joint_posterior(cluster.gwas_configuration_posteriors)

	# Organise information into a hash 
	res = dict((config, config_sample.posterior[config_sample.configurations[config]]) for config in config_sample.configurations)
	res['_CLUSTER'] = sum_posteriors

	return res


def sign(number):
	"""
			Returns the sign of the number (-1, 0 or 1)
			Arg1: float
			Returntype: int
	"""
	if number > 0:
		return 1
	elif number < 0:
		return -1
	else:
		return 0


def z_score_from_pvalue(p_value, direction):
	"""
			Estimates z-score from p-value and effect direction
			Arg1: float
			Arg2: float
			Returntype: float
	"""
	global scipy
	if scipy is None:
		import scipy
		import scipy.stats
	if p_value == 0:
		p_value = 4.2e-317
	return -scipy.stats.norm.ppf(p_value/2) * sign(direction)
