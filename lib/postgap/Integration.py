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


import collections
import math
import os
import sys
import os.path
import cPickle as pickle
from postgap.DataModel import *
import postgap.GWAS
import postgap.Cisreg
import postgap.Reg
import postgap.REST
import postgap.LD
import postgap.EFO
import postgap.Ensembl_lookup
import postgap.Globals
import postgap.RegionFilter
import postgap.FinemapIntegration
import postgap.Finemap
from postgap.Utils import *
import logging
import sys
phenotype_cache = ()


def diseases_to_genes(diseases, efos, population, tissues):
	"""

			Associates genes from a list of diseases
			Args:
			* [ string ] (trait descriptions - free strings)
			* [ string ] (trait EFO identifiers)
			* string (population name)
			* [ string ] (tissue names)
			Returntype: [ GeneCluster_Association ]

	"""
	return gwas_snps_to_genes(diseases_to_gwas_snps(diseases, efos), population, tissues)


def diseases_to_gwas_snps(diseases, efos):
	"""

			Associates gwas_snps from a list of diseases
			Args:
			* [ string ] (trait descriptions - free strings )
			* [ string ] (trait EFO identifiers)
			Returntype: [ GWAS_SNP ]

	"""
	res = filter(lambda X: X.pvalue < postgap.Globals.GWAS_PVALUE_CUTOFF,
				 scan_disease_databases(diseases, efos))

	logging.info("Found %i GWAS SNPs associated to diseases (%s) or EFO IDs (%s) after p-value filter (%f)" %
				 (len(res), ", ".join(diseases), ", ".join(efos), postgap.Globals.GWAS_PVALUE_CUTOFF))

	return res


def scan_disease_databases(diseases, efos):
	"""

			Associates gwas_snps from a list of diseases
			Args:
			* [ string ] (trait descriptions)
			* [ string ] (trait EFO identifiers)
			Returntype: [ GWAS_SNP ]

	"""
	if postgap.Globals.GWAS_SUMMARY_STATS_FILE is not None:
		gwas_associations = postgap.GWAS.GWAS_File().run(diseases, efos)
	elif postgap.Globals.GWAS_adaptors == None:
		logging.info("Searching for GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in all databases" % (
			", ".join(diseases), ", ".join(efos)))
		gwas_associations = concatenate(source().run(
			diseases, efos) for source in postgap.GWAS.sources)
	else:
		logging.info("Searching for GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in (%s)" % (
			", ".join(diseases), ", ".join(efos), ", ".join(postgap.Globals.GWAS_adaptors)))
		gwas_associations = concatenate(source().run(
			diseases, efos) for source in postgap.GWAS.get_filtered_subclasses(postgap.Globals.GWAS_adaptors))

	logging.info("Done searching for GWAS SNPs. Found %s gwas associations." % len(
		gwas_associations))

	associations_by_snp = dict()
	for gwas_association in gwas_associations:
		# Sanity filter to avoid breaking downstream code
		# Looking at you GWAS DB, "p-value = 0", pshaw!
		if gwas_association.pvalue <= 0:
			continue
		if gwas_association.sample_size <= 0:
			continue

		if postgap.Globals.PERFORM_BAYESIAN:
			if gwas_association.odds_ratio is not None:
				z_score = postgap.FinemapIntegration.z_score_from_pvalue(
					gwas_association.pvalue, float(gwas_association.odds_ratio) - 1)
				beta = math.log(gwas_association.odds_ratio)
			elif gwas_association.beta_coefficient is not None:
				z_score = postgap.FinemapIntegration.z_score_from_pvalue(
					gwas_association.pvalue, gwas_association.beta_coefficient)
				beta = gwas_association.beta_coefficient
			else:
				continue
		else:
			z_score = None
			beta = None

		if gwas_association.snp not in associations_by_snp or associations_by_snp[gwas_association.snp].pvalue < gwas_association.pvalue:
			associations_by_snp[gwas_association.snp] = GWAS_SNP(
				snp=SNP(rsID=gwas_association.snp, chrom=None,
						pos=None, approximated_zscore=None),
				pvalue=gwas_association.pvalue,
				evidence=[gwas_association],
				z_score=z_score,
				beta=beta
			)

	gwas_snps = associations_by_snp.values()

	if postgap.Globals.GWAS_adaptors == None:
		logging.info("Found %i unique GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in all databases" % (
			len(gwas_snps), ", ".join(diseases), ", ".join(efos)))
	else:
		logging.info("Found %i unique GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in (%s)" % (
			len(gwas_snps), ", ".join(diseases), ", ".join(efos), ", ".join(postgap.Globals.GWAS_adaptors)))

	return gwas_snps

def gwas_snps_to_genes(gwas_snps, population, tissue_weights):
	"""

			Associates Genes to gwas_snps of interest
			Args:
			* [ GWAS_Association ]
			* string (population name)
			* { tissue_name: scalar (weight) }
			Returntype: [ GeneCluster_Association ]

	"""
	# Must set the tissue settings before separating out the gwas_snps
	# TODO: tissue_weights may be deprecated once ML weighting is implemented
	if tissue_weights is None:
		tissue_weights = gwas_snps_to_tissue_weights(gwas_snps)

	if gwas_snps is None and postgap.Globals.CLUSTER_FILE is not None:
		cluster_fn = postgap.Globals.CLUSTER_FILE
		infile = open(cluster_fn, 'rb')
		cluster = pickle.load(infile)
		infile.close()

		# Perform GWAS finemapping of the clusters
		if postgap.Globals.PERFORM_BAYESIAN:
			logging.info("\tperform GWAS finemapping for cluster %s" % (cluster_fn))
		
			try:
				cluster = postgap.FinemapIntegration.finemap_gwas_cluster(cluster, population)
			except:
				logging.info("\tcluster %s failed finemapping" % (cluster_fn))
				sys.exit(0) # should quit here??? Reference Integration.py - L242-246

		res = cluster_to_genes(cluster, tissue_weights, population)

		logging.info("\tFound %i genes associated to cluster %s" % (len(res), cluster_fn))
	elif gwas_snps is not None:
		clusters = cluster_gwas_snps(gwas_snps, population)

		if postgap.Globals.CLUSTER_DIR is not None:
			# create pickle files in directory cluster_dir
			for cluster in clusters:
				# create unique sample label for each (raw) cluser
				chrom = cluster.ld_snps[0].chrom
				start = min(ld_snp.pos for ld_snp in cluster.ld_snps)
				end = max(ld_snp.pos for ld_snp in cluster.ld_snps)

				# create pickle files in directory cluster_dir		
				cluster_fn = 'GWAS_Cluster_%s:%i-%i' % (chrom, start, end)
				outfile = open(postgap.Globals.CLUSTER_DIR + cluster_fn, 'wb')
				pickle.dump(cluster, outfile)
				outfile.close()
			logging.info("save cluster info into a file and quit")
			sys.exit(0)

		# Perform GWAS finemapping of the clusters
		if postgap.Globals.PERFORM_BAYESIAN:
			logging.info("\tperform GWAS finemapping for all clusters")
		
			finemap_clusters = []
			for cluster in clusters:
				try:
					finemap_clusters.append(postgap.FinemapIntegration.finemap_gwas_cluster(cluster, population))
				except:
					continue
				
			if not finemap_clusters:
				logging.info("all the clusters failed finemapping")
				sys.exit(0) # should quit here??? or 'finemap_clusters = clusters'
		
			for cluster in finemap_clusters:
				for gwas_snp in cluster.gwas_snps:
					assert gwas_snp.snp.rsID in [ld_snp.rsID for ld_snp in cluster.ld_snps]
		
			clusters = finemap_clusters

		res = concatenate(cluster_to_genes(cluster, tissue_weights, population) for cluster in clusters)

		logging.info("\tFound %i genes associated to all clusters" % (len(res)))

	if len(res) == 0:
		return []
	else:
		return sorted(res, key=lambda X: X.score)

def clusters_to_genes(clusters, population, tissue_weights):
	"""

			Associates genes to a set of clusters 
			Args:
			* [ Cluster ]
			* string (population name)
			* {tissue_name: scalar (weight) }
			Returntype: [ string ]

	"""
	# Collect regulatory and cis-regulatory evidence across clusters
	cluster_associations = [(cluster, ld_snps_to_genes(
		cluster.ld_snps, tissue_weights)) for cluster in clusters]
	with open('cluster_association.pkl', 'w') as f:
		pickle.dump(cluster_associations, f)
	# with open('cluster_association.pkl') as f:
	#		cluster_associations= pickle.load(f)
	# If required, perform genome-wide GWAS finemapping

	# Extract the GWAS and eQTL F.A parameters (lambdas)
	# TODO: Remove magic filepaths
	# with open('GWAS_lambdas_'+os.path.basename(postgap.Globals.GWAS_SUMMARY_STATS_FILE), 'w') as fw1:
	#	 fw1.write('cluster\tsource\tlambdas\n')
	# with open('eQTL_lambdas_'+os.path.basename(postgap.Globals.GWAS_SUMMARY_STATS_FILE), 'w') as fw2:
	#	 fw2.write('cluster\ttissue\tgene\tlambdas\n')
	with open(postgap.Globals.OUTPUT+'_GWAS_lambdas.txt', 'w') as fw1:
		fw1.write('source\tlambdas\n')
	with open(postgap.Globals.OUTPUT+'_eQTL_lambdas.txt', 'w') as fw2:
		fw2.write('cluster\ttissue\tgene\tlambdas\n')
	# ===

	if postgap.Globals.PERFORM_BAYESIAN:
		cluster_associations = postgap.FinemapIntegration.compute_gwas_posteriors(
			cluster_associations, population)
	# print cluster_associations[0][0].ld_matrix

	# Perform cluster by cluster finemapping
	res = concatenate(cluster_to_genes(cluster, associations, tissue_weights, population)
					  for cluster, associations in cluster_associations)

	logging.info("\tFound %i genes associated to all clusters" % (len(res)))

	if len(res) == 0:
		return []
	else:
		return sorted(res, key=lambda X: X.score)


def gwas_snps_to_tissue_weights(gwas_snps):
	"""

			Associates list of tissues to list of gwas_snps
			Args:
			* [ GWAS_SNP ]
			Returntype: [ string ]

	"""
	return ['Whole_Blood']  # See FORGE??


def cluster_gwas_snps(gwas_snps, population):
	"""

			Bundle together gwas_snps within LD threshold
			* [ GWAS_SNP ]
			* String (population name)
			Returntype: [ GWAS_Cluster ]

	"""
	gwas_snp_locations = get_gwas_snp_locations(gwas_snps)

	logging.info("Found %i locations from %i GWAS SNPs" %
				 (len(gwas_snp_locations), len(gwas_snps)))

	# For every gwas snp location, create the preclusters by simple LD expansion of independent SNPs.
	#
	preclusters = filter(lambda X: X is not None, [gwas_snp_to_precluster(
		gwas_snp_location, population) for gwas_snp_location in gwas_snp_locations])

	# Remove preclusters having a SNP in a blacklisted region
	#
	filtered_preclusters = postgap.RegionFilter.region_filter(preclusters)
	for precluster in filtered_preclusters:
		for gwas_snp in precluster.gwas_snps:
			assert gwas_snp.snp.rsID in [ld_snp.rsID for ld_snp in precluster.ld_snps]
	
	# merge clusters using the new strategy and remove clusters having less than 10 ld_snps as they are less informative
	processed_preclusters = remove_overlaps(merge_preclusters_ld(filtered_preclusters))
	raw_clusters = remove_tiny_clusters(processed_preclusters)
	for cluster in raw_clusters:
		for gwas_snp in cluster.gwas_snps:
			assert gwas_snp.snp.rsID in [ld_snp.rsID for ld_snp in cluster.ld_snps]

	logging.info("Found %i clusters from %i GWAS SNP locations" % (len(raw_clusters), len(gwas_snp_locations)))

	return raw_clusters

def get_gwas_snp_locations(gwas_snps):
	"""

			Extract locations of GWAS SNP:
			Args
			* GWAS_SNP
			Returntype: [ GWAS_SNP ]

	"""
	original_gwas_snp = dict((gwas_snp.snp.rsID, gwas_snp)
							 for gwas_snp in gwas_snps)
	mapped_snps = postgap.Ensembl_lookup.get_snp_locations(
		original_gwas_snp.keys())

	return [
		GWAS_SNP(
			snp=mapped_snp,
			pvalue=original_gwas_snp[mapped_snp.rsID].pvalue,
			evidence=original_gwas_snp[mapped_snp.rsID].evidence,
			z_score=original_gwas_snp[mapped_snp.rsID].z_score,
			beta=original_gwas_snp[mapped_snp.rsID].beta
		)
		for mapped_snp in mapped_snps
		if mapped_snp.rsID in original_gwas_snp
	]

def distance_between_preclusters(A, B):
	min_distance = None
	for SNPA in A.gwas_snps:
		for SNPB in B.gwas_snps:
			if SNPA.snp.chrom == SNPB.snp.chrom and (min_distance is None or min_distance > abs(SNPA.snp.pos - SNPB.snp.pos)):
				min_distance = abs(SNPA.snp.pos - SNPB.snp.pos)
	return min_distance


def merge_preclusters_ld(preclusters):
	"""
		Bundle together preclusters that satisfy the two criteria at the same time:
		* 1. gwas_snp of clusterA is within ld_snps of clusterB (in LD)
		* 2. (p-value of gwas_snp of clusterB * 10-3) <= p-value of gwas_snp of clusterA (<= p-value of gwas_snp of clusterB)
		Args:
		* [ Cluster ]
		Returntype: [ Cluster ]
	"""
	# sort preclusters by the order of p-value of GWAS SNPs, from lowest to highest
	preclusters.sort(key=lambda cluster: cluster.gwas_snps[0].pvalue)
	
	# create a vaiable for the status of each precluster whether it has been grouped with any other
	group_status = [False] * len(preclusters)
	
	merged_clusters = []
	
	# clusterA - take its gwas_snp for checking
	for i,clusterA in enumerate(preclusters[ :-1]):
		# generate temp_cluster
		if group_status[i]:
			where_is_A = [n for n,cluster in enumerate(merged_clusters) if clusterA.gwas_snps[0].snp.rsID in [gwas_snp.snp.rsID for gwas_snp in cluster.gwas_snps]]
			assert len(where_is_A), 'clusterA should be in one place only in merged_clusters'
			
			temp_cluster = merged_clusters.pop(where_is_A[0])
		else:
			temp_cluster = None
		
		# get the index of what to be merged with clusterA
		what_to_merge = []
		
		# clusterB - take its ld_snps and p-value of gwas_snp for comparison
		for k,clusterB in enumerate(preclusters[i + 1: ]):
			k = k + i + 1
			
			# two have the posibility to be merged only when they are from the same chromosome
			if clusterB.ld_snps[0].chrom == clusterA.ld_snps[0].chrom:
				# first, check p-value
				if clusterB.gwas_snps[0].pvalue * 1e-3 <= clusterA.gwas_snps[0].pvalue:
					# secondly, check ld_snps
					if clusterA.gwas_snps[0].snp.rsID in [ld_snp.rsID for ld_snp in clusterB.ld_snps]:
						# thirdly, check whether in temp_cluster
						if temp_cluster is None:
							what_to_merge += [k]
						else:
							if clusterB.gwas_snps[0].snp.rsID not in [gwas_snp.snp.rsID for gwas_snp in temp_cluster.gwas_snps]:
								what_to_merge += [k]
				# no more following clusters will meet the criterion, so break out of clusterB's for-loop
				else:
					break
		
		# if nothing to merge, put temp_cluster back
		if len(what_to_merge) == 0:
			if temp_cluster is not None:
				merged_clusters.append(temp_cluster)
		# if something to merge, do merge_clusters
		else:
			clusters_to_merge = []
			
			# for any has been merged with clusters other than clusterA, get that merged one and merge it with clusterA
			if any(group_status[m] for m in what_to_merge):
				for n in set(n for m in what_to_merge if group_status[m] for n,cluster in enumerate(merged_clusters) if preclusters[m].gwas_snps[0].snp.rsID in [gwas_snp.snp.rsID for gwas_snp in cluster.gwas_snps]):
					clusters_to_merge.append(merged_clusters.pop(n))
			
			# for any has not been merged with others, get it to merge with clusterA
			for m in what_to_merge:
				if not group_status[m]:
					clusters_to_merge.append(preclusters[m])
			
			temp_cluster = merge_clusters([clusterA if temp_cluster is None else temp_cluster] + clusters_to_merge)
			merged_clusters.append(temp_cluster)
			
			# update the status of these clusters after merging
			for m in [i] + what_to_merge:
				group_status[m] = True
	
	# add those ungrouped clusters into merged_clusters
	clusters = merged_clusters + [preclusters[n] for n,s in enumerate(group_status) if not s]
	
	return clusters

def merge_clusters(cluster_list):
	"""
		Merges a list of clusters into a single one:
		Args:
		* [ Cluster ], eg: [clusterA] + [clusterB]
		Returntype: Cluster
	"""
	# Create a unique set of gwas_snps from elements in cluster_list
	merged_gwas_snps = dict((gwas_snp.snp.rsID, gwas_snp) for cluster in cluster_list for gwas_snp in cluster.gwas_snps).values()
	
	# Create a unique set of ld_snps from elements in cluster_list
	merged_ld_snps = dict((ld_snp.rsID, ld_snp) for cluster in cluster_list for ld_snp in cluster.ld_snps).values()
	
	# Create a new Cluster
	return GWAS_Cluster(
		gwas_snps=merged_gwas_snps,
		ld_snps=merged_ld_snps,
		ld_matrix=None,
		z_scores=None,
		betas=None,
		mafs=None,
		annotations=None,
		gwas_configuration_posteriors=None
	)

def remove_overlaps(preclusters):
	"""
		Remove preclusters that are less significant p-value in the two clusters who:
		* 1. share any SNPs
		* 2. have full or partial overlaps
		Args:
		* [ Cluster ]
		Returntype: [ Cluster ]
	"""
	# create a variable for the status of each cluster whether it would be kept
	keep_status = [True] * len(preclusters)
	
	# there could be some complicated situations that two clusters have the same minimum p-value, and that would call for further operations
	# so first create a variable to record those clusters
	complist = []
	
	for i,clusterA in enumerate(preclusters):
		for k,clusterB in enumerate(preclusters[i + 1: ]):
			k = k + i + 1
			
			# if clusterB is on the same chromosome as clusterA, check for overlapping
			if clusterB.ld_snps[0].chrom == clusterA.ld_snps[0].chrom:
				# get a list of all the SNPs in clusterA and a list for clusterB
				snplistA = [gwas_snp.snp.rsID for gwas_snp in clusterA.gwas_snps] + [ld_snp.rsID for ld_snp in clusterA.ld_snps]
				snplistB = [gwas_snp.snp.rsID for gwas_snp in clusterB.gwas_snps] + [ld_snp.rsID for ld_snp in clusterB.ld_snps]
				
				# if they share any SNP (in LD with each other) or share any overlaps on the coordinates
				if any(snp in snplistA for snp in snplistB) or not (max(ld_snp.pos for ld_snp in clusterB.ld_snps) < min(ld_snp.pos for ld_snp in clusterA.ld_snps)  or min(ld_snp.pos for ld_snp in clusterB.ld_snps) > max(ld_snp.pos for ld_snp in clusterA.ld_snps)):
					# keep the cluster with more significant p-value
					if min(gwas_snp.pvalue for gwas_snp in clusterA.gwas_snps) < min(gwas_snp.pvalue for gwas_snp in clusterB.gwas_snps):
						keep_status[k] = False
					elif min(gwas_snp.pvalue for gwas_snp in clusterA.gwas_snps) > min(gwas_snp.pvalue for gwas_snp in clusterB.gwas_snps):
						keep_status[i] = False
					else:
						if keep_status[i] and keep_status[k]:
							complist.append((i, k))
	
	# do something about those complicated clusters, if there are any
	#	1. if they share a SNP and have the same minimum p-value, they should be merged.
	#	2. if they do not share a SNP, but overlap, we leave them as they are.
	merged_clusters = []

	if len(complist) > 0 and any(keep_status[c] for t in complist for c in t):
		for t in complist:
			# if both clusters have NOT been rejected by any clusters
			if keep_status[t[0]] and keep_status[t[1]]:
				assert min(gwas_snp.pvalue for gwas_snp in clusterA.gwas_snps) == min(gwas_snp.pvalue for gwas_snp in clusterB.gwas_snps), 'this is a bug!'
				
				clusterA =  preclusters[t[0]]
				clusterB =  preclusters[t[1]]
				snplistA = [gwas_snp.snp.rsID for gwas_snp in clusterA.gwas_snps] + [ld_snp.rsID for ld_snp in clusterA.ld_snps]
				snplistB = [gwas_snp.snp.rsID for gwas_snp in clusterB.gwas_snps] + [ld_snp.rsID for ld_snp in clusterB.ld_snps]
				
				# merge two clusters if they any share SNP(s) and have the same min(p-value)
				if any(snp in snplistA for snp in snplistB):
					merged_clusters.append(merge_clusters([clusterA, clusterB]))
	
					keep_status[t[0]] = False
					keep_status[t[1]] = False
	
	# clear out removed clusters, and append merged clusters
	clusters = [cluster for s, cluster in zip(keep_status, preclusters) if s] + merged_clusters
	
	return clusters

def remove_tiny_clusters(preclusters):
	"""
		Remove preclusters that have less than 10 SNPs in them
		* [ Cluster ]
		Returntype: [ Cluster ]
	"""
	clusters = [cluster for cluster in preclusters if len([ld_snp.rsID for ld_snp in cluster.ld_snps]) >= 10]
	
	return clusters

def cluster_to_genes(cluster, tissues, population):
	"""
		Associated Genes to a cluster of gwas_snps
		Args:
		* [ Cluster ]
		* { tissue_name: scalar (weights) }
		* string (population name)
		Returntype: [ GeneCluster_Association ]
	"""
	if postgap.Globals.PERFORM_BAYESIAN:
		assert len(cluster.ld_snps) == cluster.ld_matrix.shape[0], (len(cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
		assert len(cluster.ld_snps) == cluster.ld_matrix.shape[1], (len(cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
		# Obtain interaction data from LD snps
		associations, eQTL_hash = ld_snps_to_genes([snp for snp in cluster.ld_snps], tissues)
		gene_tissue_posteriors = postgap.FinemapIntegration.compute_joint_posterior(cluster, eQTL_hash)

		res = [
			GeneCluster_Association(
				gene = gene,
				score = None,
				collocation_posterior = gene_tissue_posteriors[gene],
				cluster = cluster,
				evidence = filter(lambda X: X.gene == gene, associations), # This is a [ GeneSNP_Association ]
				eQTL_hash = eQTL_hash,
				r2 = None
			)
			for gene in gene_tissue_posteriors
		]

		logging.info("\tFound %i genes associated" % (len(res)))

	else:
		# Compute LD from top SNP
		top_gwas_hit = sorted(cluster.gwas_snps, key=lambda X: X.pvalue)[-1]
		ld = postgap.LD.get_lds_from_top_gwas(top_gwas_hit.snp, cluster.ld_snps, population=population)

		# Obtain interaction data from LD snps
		associations, eQTL_hash = ld_snps_to_genes([snp for snp in cluster.ld_snps if snp in ld], tissues)
		
		# Compute gene score
		gene_scores = dict(
			((association.gene, association.snp), (association, association.score * ld[association.snp]))
			for association in associations
		)

		if len(gene_scores) == 0: return []

		# OMIM exception
		max_score = max(X[1] for X in gene_scores.values())
		for gene, snp in gene_scores:
			if len(gene_to_phenotypes(gene)):
				gene_scores[(gene, snp)][1] = max_score

		# PICS score precomputed and normalised
		pics = PICS(ld, top_gwas_hit.pvalue)

		# Compute posterior
		res = [
			GeneCluster_Association(
				gene = gene,
				score = total_score(pics[snp], gene_scores[(gene, snp)][1]),
				collocation_posterior = None,
				cluster = cluster,
				evidence = gene_scores[(gene, snp)][:1], # This is a [ GeneSNP_Association ]
				eQTL_hash = eQTL_hash,
				r2 = ld[snp]
			)
			for (gene, snp) in gene_scores
			if snp in pics
		]

		logging.info("\tFound %i genes associated around GWAS SNP %s" % (len(res), top_gwas_hit.snp.rsID))

	# Pick the association with the highest score
	return sorted(res, key=lambda X: X.score)

def PICS(ld, pvalue):
	minus_log_pvalue = - math.log(pvalue) / math.log(10)
	SD = dict()
	Mean = dict()
	prob = dict()
	sum = 0

	for snp in ld.keys():
		if snp in ld:
			SD[snp] = math.sqrt(1 - math.sqrt(ld[snp]) ** 6.4) * \
				math.sqrt(minus_log_pvalue) / 2
			Mean[snp] = ld[snp] * minus_log_pvalue
		else:
			# Defaults for remote SNPs
			SD[snp] = 0
			Mean[snp] = 1 + minus_log_pvalue

		# Calculate the probability of each SNP
		if SD[snp]:
			prob[snp] = 1 - pnorm(minus_log_pvalue, Mean[snp], SD[snp])
		else:
			prob[snp] = 1

		# Normalisation sum
		sum += prob[snp]

	# Normalize the probabilies so that their sum is 1.
	return dict((snp, prob[snp] / sum) for snp in prob.keys())


def pnorm(x, mu, sd):
	"""

			Normal distribution PDF
			Args:
			* scalar: variable
			* scalar: mean
			* scalar: standard deviation
			Return type: scalar (probability density)

	"""
	return math.exp(- ((x - mu) / sd) ** 2 / 2) / (sd * 2.5)


def total_score(pics, gene_score):
	"""

			Computes a weird mean function from ld_snp PICs score and Gene/SNP association score
			Args:
			* PICS: scalar
			* gene_score: scalar
			Returntype: scalar

	"""
	if pics is None:
		return None
	A = pics * (pics ** (1/3))
	B = gene_score * (gene_score ** (1/3))
	return ((A + B) / 2) ** 3


def rsIDs_to_genes(snp, tissues):
	"""

			Associates genes to single SNP
			Args:
			* SNP
			* [ string ] (tissues)
			Returntype: [ GeneSNP_Association ]

	"""
	if tissues is None:
		tissues = ["Whole_Blood"]
	return ld_snps_to_genes(postgap.Ensembl_lookup.get_snp_locations([snp]), tissues)


def ld_snps_to_genes(ld_snps, tissues):
	"""
		Associates genes to LD linked SNPs
		Args:
		* [ SNP ]
		* [ string ] (tissues)
		* Dict SNP => float
		Returntype:
		* [ GeneCluster_Association ]
		* Hash of hashes: Gene => Tissue => SNP => Float
	"""
	# Search for SNP-Gene pairs:
	cisreg = cisregulatory_evidence(ld_snps, tissues) # Hash of hashes: SNP => Gene => Cisregulatory_Evidence

	# Organise eQTL data into an easily read hash of hashes
	gene_tissue_snp_eQTL_hash = organise_eQTL_data(cisreg)
	
	# filter eQTL data that goes into the eQTL hash
	cisreg = filter_eQTL_GTEx(cisreg)

	# Extract SNP specific info:
	reg = regulatory_evidence(cisreg.keys(), tissues) # Hash: SNP => [ Regulatory_evidence ]
		
	associations = concatenate((create_SNP_GeneSNP_Associations(snp, reg[snp], cisreg[snp]) for snp in cisreg))
	
	return associations, gene_tissue_snp_eQTL_hash

def organise_eQTL_data(cisreg):
	"""
		Organise unsorted eQTL data into easily read hash of hashes:
		Arg1: Hash of hashes: SNP => Gene => Cisregulatory_Evidence
		Returntype: Hash of hashes: Gene => Tissue => SNP => Float
	"""
	res = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(float)))
	for snp in cisreg:
		for gene in cisreg[snp]:
			for evidence in cisreg[snp][gene]:
				# save GTEx eQTL data - z-score, p-value and beta
				if evidence.source == 'GTEx':
					res[gene][evidence.tissue][snp.rsID] = (evidence.z_score, evidence.pvalue, evidence.beta)
	return res

def filter_eQTL_GTEx(cisreg):
	"""
		Filter out eQTL data that goes into the eQTL hash
	"""
	for snp in cisreg:
		for gene in cisreg[snp]:
			cisreg[snp][gene] = filter(lambda X: X.source != 'GTEx', cisreg[snp][gene])
	
	return cisreg

def create_SNP_GeneSNP_Associations(snp, reg, cisreg):
	"""

			Associates gene to LD linked SNP
			Args:
			* SNP
			* [ Regulatory_Evidence ]
			* [ Cisregulatory_Evidence ]
			Returntype: [GeneSNP_Association]

	"""
	intermediary_scores, gene_scores = compute_v2g_scores(reg, cisreg)
	rank = dict((score, index) for index, score in enumerate(
		sorted(gene_scores.values(), reverse=True)))

	return [GeneSNP_Association(
		gene=gene,
		snp=snp,
		cisregulatory_evidence=cisreg[gene],
		regulatory_evidence=reg,
		intermediary_scores=intermediary_scores[gene],
		score=gene_scores[gene],
		rank=rank[gene_scores[gene]] + 1)
		for gene in cisreg]


def compute_v2g_scores(reg, cisreg):
	"""

			Goes through evidence and scores associations to a SNP
			Args:
			* [ Regulatory_Evidence ]
			* [ Cisregulatory_Evidence ] 
			Returntype: dict(Gene: dict(string: float)), dict(Gene: float)

	"""

	intermediary_scores = dict()
	gene_scores = dict()
	for gene in cisreg:
		intermediary_scores[gene] = collections.defaultdict(int)
		seen = set()
		for evidence in cisreg[gene] + reg:
			if evidence.source not in seen or float(evidence.score) > intermediary_scores[gene][evidence.source]:
				intermediary_scores[gene][evidence.source] = float(
					evidence.score)
				seen.add(evidence.source)

			# VEP stats
			if evidence.source == 'VEP':
				intermediary_scores[gene]['VEP_count'] += 1
				intermediary_scores[gene]['VEP_sum'] += float(evidence.score)

			if evidence.source == 'GTEx':
				intermediary_scores[gene][evidence.tissue] = float(
					evidence.score)

			# Ad hoc bounds defined here:
		# PCHiC
		intermediary_scores[gene]['PCHiC'] = min(
			intermediary_scores[gene]['PCHiC'], 1)

		# VEP
		if 'VEP' in intermediary_scores[gene]:
			intermediary_scores[gene]['VEP_mean'] = intermediary_scores[gene]['VEP_sum'] / \
				intermediary_scores[gene]['VEP_count']

		# Weighted sum
		gene_scores[gene] = sum(intermediary_scores[gene][source] * postgap.Globals.EVIDENCE_WEIGHTS[source]
								for source in intermediary_scores[gene] if source in postgap.Globals.EVIDENCE_WEIGHTS)

	return intermediary_scores, gene_scores


def cisregulatory_evidence(ld_snps, tissues):
	"""

			Associates genes to LD linked SNP
			Args:
			* [ SNP ]
			* [ string ] (tissues)
			Returntype: Hash of hashes SNP => Gene => Cisregulatory_Evidence

	"""
	if postgap.Globals.Cisreg_adaptors == None:
		logging.info(
			"Searching for cis-regulatory data on %i SNPs in all databases" % (len(ld_snps)))
		evidence = concatenate(source().run(ld_snps, tissues)
							   for source in postgap.Cisreg.sources)
	else:
		logging.info("Searching for cis-regulatory data on %i SNPs in (%s)" %
					 (len(ld_snps), ", ".join(postgap.Globals.Cisreg_adaptors)))
		evidence = concatenate(source().run(ld_snps, tissues)
							   for source in postgap.Cisreg.get_filtered_subclasses(postgap.Globals.Cisreg_adaptors))

	flaten_evidence = []
	for evi in evidence:
		if type(evi) == list:
			for e in evi:
				flaten_evidence.append(e)
		else:
			flaten_evidence.append(evi)

	filtered_evidence = filter(
		lambda association: association is not None and association.gene is not None and association.gene.biotype == "protein_coding", flaten_evidence)

	# Group by snp, then gene:
	res = collections.defaultdict(lambda: collections.defaultdict(list))
	for association in filtered_evidence:
		assert type(
			association) is postgap.DataModel.Cisregulatory_Evidence, "association is Cisregulatory_Evidence"
		res[association.snp][association.gene].append(association)

	if postgap.Globals.Cisreg_adaptors == None:
		logging.debug(
			("Found %i cis-regulatory interactions in all databases" % (len(res))))
	else:
		logging.debug(("Found %i cis-regulatory interactions in (%s)" %
					   (len(res), ", ".join(postgap.Globals.Cisreg_adaptors))))
	return res


def regulatory_evidence(snps, tissues):
	"""

			Extract regulatory evidence linked to SNPs and stores them in a hash
			* [ SNP ]
			* [ string ]
			Returntype: Hash: SNP => [ Regulatory_evidence ]

	"""
	if postgap.Globals.Reg_adaptors == None:
		logging.info(
			"Searching for regulatory data on %i SNPs in all databases" % (len(snps)))
		res = concatenate(source().run(snps, tissues)
						  for source in postgap.Reg.sources)
		logging.info("Found %i regulatory SNPs among %i in all databases" % (
			len(res), len(snps)))
	else:
		logging.info("Searching for regulatory data on %i SNPs in (%s)" %
					 (len(snps), ", ".join(postgap.Globals.Reg_adaptors)))
		res = concatenate(source().run(snps, tissues)
						  for source in postgap.Reg.get_filtered_subclasses(postgap.Globals.Reg_adaptors))
		logging.info("Found %i regulatory SNPs among %i in (%s)" % (
			len(res), len(snps), ", ".join(postgap.Globals.Reg_adaptors)))

	# Group by SNP
	hash = collections.defaultdict(list)
	for hit in res:
		hash[hit.snp].append(hit)

	return hash


def gene_to_phenotypes(gene):
	"""

			Look up phenotype annotations for gene
			Args:
			* ENSG stable ID
			Return type: [ OMIM Phenotype ]

	"""
	return []  # TODO and remove stopper
	if gene not in phenotype_cache:
		phenotype_cache[gene['stable_id']
						] = postgap.Ensembl_lookup.get_gene_phenotypes(gene)

	return phenotype_cache[gene['stable_id']]


def get_all_tissues():
	server = 'http://rest.ensembl.org'
	ext = '/eqtl/tissue/homo_sapiens?content-type=application/json'
	try:
		return postgap.REST.get(server, ext)

	except:
		return None
