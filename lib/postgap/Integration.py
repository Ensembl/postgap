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
import scipy


import sys

phenotype_cache = ()

def diseases_to_genes(diseases, efos, populations, tissues):
	"""

		Associates genes from a list of diseases
		Args:
		* [ string ] (trait descriptions - free strings)
		* [ string ] (trait EFO identifiers)
		* [ string ] (population names)
		* [ string ] (tissue names)
		Returntype: [ GeneCluster_Association ]

	"""
	return gwas_snps_to_genes(diseases_to_gwas_snps(diseases, efos), populations, tissues)

def diseases_to_gwas_snps(diseases, efos):
	"""

		Associates gwas_snps from a list of diseases
		Args:
		* [ string ] (trait descriptions - free strings )
		* [ string ] (trait EFO identifiers)
		Returntype: [ GWAS_SNP ]

	"""
	res = filter(lambda X: X.pvalue < postgap.Globals.GWAS_PVALUE_CUTOFF, scan_disease_databases(diseases, efos))

	logging.info("Found %i GWAS SNPs associated to diseases (%s) or EFO IDs (%s) after p-value filter (%f)" % (len(res), ", ".join(diseases), ", ".join(efos), postgap.Globals.GWAS_PVALUE_CUTOFF))

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
		logging.info("Searching for GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in all databases" % (", ".join(diseases), ", ".join(efos)))
		gwas_associations = concatenate(source().run(diseases, efos) for source in postgap.GWAS.sources)
	else:
		logging.info("Searching for GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in (%s)" % (", ".join(diseases), ", ".join(efos), ", ".join(postgap.Globals.GWAS_adaptors)))
		gwas_associations = concatenate(source().run(diseases, efos) for source in postgap.GWAS.get_filtered_subclasses(postgap.Globals.GWAS_adaptors))

	logging.info("Done searching for GWAS SNPs. Found %s gwas associations." % len(gwas_associations) )

	# DEBUG
	for gwas_association in gwas_associations:
		snp = gwas_association.snp
		print "\t".join(map(str, [snp.chrom, snp.pos, snp.rsID, gwas_association.pvalue]))

	associations_by_snp = dict()
	for gwas_association in gwas_associations:
		# Sanity filter to avoid breaking downstream code
		# Looking at you GWAS DB, "p-value = 0", pshaw!
		if gwas_association.pvalue <= 0:
			continue
		
		if postgap.Globals.PERFORM_BAYESIAN:
			if gwas_association.odds_ratio is not None:
				z_score = postgap.FinemapIntegration.z_score_from_pvalue(gwas_association.pvalue, gwas_association.odds_ratio)
			elif gwas_association.beta_coefficient is not None:
				z_score = postgap.FinemapIntegration.z_score_from_pvalue(gwas_association.pvalue, gwas_association.beta_coefficient)
			else:
				continue
		else:
			z_score = None

		if gwas_association.snp not in associations_by_snp or associations_by_snp[gwas_association.snp].pvalue < gwas_association.pvalue:
			associations_by_snp[gwas_association.snp] = GWAS_SNP(
				snp = gwas_association.snp,
				pvalue = gwas_association.pvalue,
				evidence = [ gwas_association ],
				z_score   = z_score,
			)

	gwas_snps = associations_by_snp.values()

	if postgap.Globals.GWAS_adaptors == None:
		logging.info("Found %i unique GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in all databases" % (len(res), ", ".join(diseases), ", ".join(efos)))
	else:
		logging.info("Found %i unique GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in (%s)" % (len(res), ", ".join(diseases), ", ".join(efos), ", ".join(postgap.Globals.GWAS_adaptors)))

	return gwas_snps

def gwas_snps_to_genes(gwas_snps, populations, tissue_weights):
	"""

		Associates Genes to gwas_snps of interest
		Args:
		* [ GWAS_Association ]
		* { population_name: scalar (weight) }
		* { tissue_name: scalar (weight) }
		Returntype: [ GeneCluster_Association ]

	"""
	# Must set the tissue settings before separating out the gwas_snps
	if tissue_weights is None:
		tissue_weights = gwas_snps_to_tissue_weights(gwas_snps)

	clusters = cluster_gwas_snps(gwas_snps, populations)
	
	if len(clusters)>0:
		
		from postgap.Globals import finemap_gwas_clusters_directory
		logging.info("Writing %i clusters to %s" % (len(clusters), finemap_gwas_clusters_directory))
		write_gwas_clusters_to_files(clusters, finemap_gwas_clusters_directory)

	return clusters
	

def clusters_to_genes(clusters, populations, tissue_weights):
	"""

		Associates genes to a set of clusters 
		Args:
		* [ Cluster ]
		* { population_name: scalar (weight) }
		* {tissue_name: scalar (weight) }
		Returntype: [ string ]

	"""
	res = concatenate(cluster_to_genes(cluster, tissue_weights, populations) for cluster in clusters)

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
	return ['Whole_Blood'] # See FORGE??

def cluster_gwas_snps(gwas_snps, populations):
	"""

		Bundle together gwas_snps within LD threshold
		* [ GWAS_SNP ]
		* [ Bio::Ensembl::Variation::Population ]
		Returntype: [ GWAS_Cluster ]

	"""
	gwas_snp_locations = get_gwas_snp_locations(gwas_snps)

	logging.info("Found %i locations from %i GWAS SNPs" % (len(gwas_snp_locations), len(gwas_snps)))

	# For every gwas snp location, create the preclusters by simple LD expansion of independent SNPs.
	#
	preclusters = filter (lambda X: X is not None, [ gwas_snp_to_precluster(gwas_snp_location, populations) for gwas_snp_location in gwas_snp_locations ])
	for precluster in preclusters:
		for gwas_snp in precluster.gwas_snps:
			assert gwas_snp.snp.rsID in [ld_snp.rsID for ld_snp in precluster.ld_snps]
	
	# Remove preclusters having a SNP in a blacklisted region
	#
	filtered_preclusters = postgap.RegionFilter.region_filter(preclusters)
	for precluster in filtered_preclusters:
		for gwas_snp in precluster.gwas_snps:
			assert gwas_snp.snp.rsID in [ld_snp.rsID for ld_snp in precluster.ld_snps]
	
	# Merge precluster that share one or more GWAS_Cluster.ld_snps.
	#
	raw_clusters = merge_preclusters_ld(merge_preclusters_distance(filtered_preclusters))
	for cluster in filtered_preclusters:
		for gwas_snp in cluster.gwas_snps:
			assert gwas_snp.snp.rsID in [ld_snp.rsID for ld_snp in cluster.ld_snps]

	if postgap.Globals.PERFORM_BAYESIAN:
		# Perform GWAS finemapping of the clusters
		# 
		clusters = [postgap.FinemapIntegration.finemap_gwas_cluster(cluster, populations) for cluster in raw_clusters]
	else:
		clusters = raw_clusters
	
	logging.info("Found %i clusters from %i GWAS SNP locations" % (len(clusters), len(gwas_snp_locations)))

	for cluster in filtered_preclusters:
		for gwas_snp in cluster.gwas_snps:
			assert gwas_snp.snp.rsID in [ld_snp.rsID for ld_snp in cluster.ld_snps]

	return clusters

def gwas_snp_to_precluster(gwas_snp, populations):
    """

        Extract neighbourhood of GWAS snp
        Args:
        * [ GWAS_SNP ]
        * [ string ] (populations)
        Returntype: GWAS_Cluster

    """
    mapped_ld_snps = postgap.LD.calculate_window(gwas_snp.snp)
    logging.info("Found %i SNPs in the vicinity of %s" % (len(mapped_ld_snps), gwas_snp.snp.rsID))
    return GWAS_Cluster(
		gwas_snps = [ gwas_snp ],
		ld_snps = mapped_ld_snps,
		ld_matrix = None,
		z_scores = None,
		gwas_configuration_posteriors = None
	)

def get_gwas_snp_locations(gwas_snps):
	"""

		Extract locations of GWAS SNP:
		Args
		* GWAS_SNP
		Returntype: [ GWAS_SNP ]

	"""
	original_gwas_snp = dict((gwas_snp.snp.rsID, gwas_snp) for gwas_snp in gwas_snps)
	mapped_snps = postgap.Ensembl_lookup.get_snp_locations(original_gwas_snp.keys())
	
	return [
		GWAS_SNP(
			snp = mapped_snp,
			pvalue = original_gwas_snp[mapped_snp.rsID].pvalue,
			evidence = original_gwas_snp[mapped_snp.rsID].evidence,
			z_score   = original_gwas_snp[mapped_snp.rsID].z_score,
		)
		for mapped_snp in mapped_snps
		if mapped_snp.rsID in original_gwas_snp
	]

def merge_preclusters_distance(preclusters):
	"""

		Bundle together preclusters that are within 100Kb
		* [ Cluster ]
		Returntype: [ Cluster ]

	"""
	input = list(preclusters)
	output = []

	for new_precluster in input:
		merged = False
		for cluster in output:
			distance = distance_between_preclusters(cluster, new_precluster)
			if distance is not None and distance < 1e5:
				output.remove(cluster)
				output.append(merge_clusters(cluster, new_precluster))
				merged = True
				break
		if not merged:
			output.append(new_precluster)

	for cluster in output:
		for gwas_snp in cluster.gwas_snps:
			assert gwas_snp.snp.rsID in [ld_snp.rsID for ld_snp in cluster.ld_snps]
	return output

def distance_between_preclusters(A, B):
	min_distance = None
	for SNPA in A.gwas_snps:
		for SNPB in B.gwas_snps:
			if SNPA.snp.chrom == SNPB.snp.chrom and (min_distance is None or min_distance > abs(SNPA.snp.pos - SNPB.snp.pos)):
				min_distance = abs(SNPA.snp.pos - SNPB.snp.pos)
	return min_distance

def merge_preclusters_ld(preclusters):
	"""

		Bundle together preclusters that share one LD snp
		* [ Cluster ]
		Returntype: [ Cluster ]

	"""
	clusters = list(preclusters)
	
	for cluster in clusters:
		chrom = cluster.gwas_snps[0].snp.chrom
		start = min(gwas_snp.snp.pos for gwas_snp in cluster.gwas_snps)
		end = max(gwas_snp.snp.pos for gwas_snp in cluster.gwas_snps)
		print "%s\t%i\t%i\t%i" % (chrom, start, end, len(cluster.gwas_snps))
	print '>>>>>>>>>>>>>>>>>>>>'
	# A dictionary that maps from snp to merged clusters
	snp_owner = dict()
	for cluster in preclusters:
		for ld_snp in cluster.ld_snps:
			# If this SNP has been seen in a different cluster
			if ld_snp in snp_owner and snp_owner[ld_snp] is not cluster:
				# Set other_cluster to that different cluster
				other_cluster = snp_owner[ld_snp]

				merged_cluster = merge_clusters(cluster, other_cluster)
				
				#    Remove the two previous clusters and replace them with 
				#    the merged cluster
				clusters.remove(cluster)
				clusters.remove(other_cluster)
				clusters.append(merged_cluster)

				#    Set the new cluster as the owner of these SNPs.
				for snp in merged_cluster.ld_snps:
					snp_owner[snp] = merged_cluster 
				for snp in cluster.ld_snps:
					snp_owner[snp] = merged_cluster 

				# Skip the rest of this cluster.
				break
			else:
				snp_owner[ld_snp] = cluster

	for cluster in clusters:
		chrom = cluster.gwas_snps[0].snp.chrom
		start = min(gwas_snp.snp.pos for gwas_snp in cluster.gwas_snps)
		end = max(gwas_snp.snp.pos for gwas_snp in cluster.gwas_snps)
		print "%s\t%i\t%i\t%i" % (chrom, start, end, len(cluster.gwas_snps))

	logging.info("\tFound %i clusters from the GWAS peaks" % (len(clusters)))

	return clusters

def merge_clusters(cluster, other_cluster):
	"""
		Merges two clusters into a single one:
		Arg1: Cluster
		Arg2: Cluster
		Returntype: Cluster
	"""
	# Merge data from current cluster into previous cluster
	merged_gwas_snps = other_cluster.gwas_snps + cluster.gwas_snps
	
	#   Create a unique set of ld snps from the current cluster 
	#   and the other cluster.
	merged_ld_snps = dict((ld_snp.rsID, ld_snp) for ld_snp in cluster.ld_snps + other_cluster.ld_snps).values()
	
	#   Create a new Cluster
	return GWAS_Cluster(
		gwas_snps = merged_gwas_snps,
		ld_snps = merged_ld_snps,
		ld_matrix = None,
		z_scores = None,
		gwas_configuration_posteriors = None
	)

def cluster_to_genes(cluster, tissues, populations):
    """

        Associated Genes to a cluster of gwas_snps
        Args:
        * [ Cluster ]
        * { tissue_name: scalar (weights) }
        * { population_name: scalar (weight) }
        Returntype: [ GeneCluster_Association ]

    """
    assert len(cluster.ld_snps) == cluster.ld_matrix.shape[0], (len(cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
    assert len(cluster.ld_snps) == cluster.ld_matrix.shape[1], (len(cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])

    if postgap.Globals.PERFORM_BAYESIAN:
      # Obtain interaction data from LD snps
      associations = ld_snps_to_genes([snp for snp in cluster.ld_snps], tissues)
      gene_tissue_posteriors = postgap.FinemapIntegration.compute_joint_posterior(cluster, associations)

      res = [
	GeneCluster_Association(
	    gene = gene,
	    score = None,
	    collocation_posterior = gene_tissue_posteriors[gene],
	    cluster = cluster,
	    evidence = filter(lambda X: X.gene == gene, associations), # This is a [ GeneSNP_Association ]
	    r2 = None
	)
	for gene in gene_tissue_posteriors
      ]

      logging.info("\tFound %i genes associated" % (len(res)))

    else:
      # Compute LD from top SNP
      top_gwas_hit = sorted(cluster.gwas_snps, key=lambda X: X.pvalue)[-1]
      ld = postgap.LD.get_lds_from_top_gwas(top_gwas_hit.snp, cluster.ld_snps, populations)

      # Obtain interaction data from LD snps
      associations = ld_snps_to_genes([snp for snp in cluster.ld_snps if snp in ld], tissues)
      
      # Compute gene score
      gene_scores = dict(
	  ((association.gene, association.snp), (association, association.score * ld[association.snp]))
	  for association in associations
      )

      if len(gene_scores) == 0:
	  return []

      # OMIM exception
      max_score = max(X[1] for X in gene_scores.values())
      for gene, snp in gene_scores:
	  if len(gene_to_phenotypes(gene)):
	      gene_scores[(gene, snp)][1] = max_score

      # PICS score precomputed and normalised
      pics = PICS(ld, top_gwas_hit.pvalue)

      # Compute posterior
      assert len(cluster.ld_snps) == cluster.ld_matrix.shape[0], (len(cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
      assert len(cluster.ld_snps) == cluster.ld_matrix.shape[1], (len(cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])

      gene_tissue_posteriors = None

      res = []
      for (gene, snp) in gene_scores:
	if snp in pics:
	  if gene_tissue_posteriors is None or gene not in gene_tissue_posteriors:
	    tissue_posteriors = None
	  else:
	    tissue_posteriors = gene_tissue_posteriors[gene]

	  res.append(GeneCluster_Association(
	      gene = gene,
	      score = total_score(pics[snp], gene_scores[(gene, snp)][1]),
	      collocation_posterior = tissue_posteriors,
	      cluster = cluster,
	      evidence = gene_scores[(gene, snp)][:1], # This is a [ GeneSNP_Association ]
	      r2 = ld[snp]
	  ))

      logging.info("\tFound %i genes associated around GWAS SNP %s" % (len(res), top_gwas_hit.snp.rsID))

    # Pick the association with the highest score
    return sorted(res, key=lambda X: X.score)

def PICS(ld, pvalue):
    minus_log_pvalue = - math.log(pvalue) / math.log(10);
    SD = dict()
    Mean = dict()
    prob = dict()
    sum = 0

    for snp in ld.keys():
        if snp in ld:
            SD[snp] = math.sqrt(1 - math.sqrt(ld[snp]) ** 6.4) * math.sqrt(minus_log_pvalue) / 2;
            Mean[snp] = ld[snp] * minus_log_pvalue;
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
	return math.exp ( - ((x - mu) / sd) ** 2 / 2 ) / (sd * 2.5)

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
		Returntype: [ GeneSNP_Association ]

	"""
	# Search for SNP-Gene pairs:
	cisreg = cisregulatory_evidence(ld_snps, tissues) # Hash of hashes SNP => Gene => Cisregulatory_Evidence
	
	# Extract SNP specific info:
	reg = regulatory_evidence(cisreg.keys(), tissues) # Hash: SNP => [ Regulatory_evidence ]
		
	return concatenate((create_SNP_GeneSNP_Associations(snp, reg[snp], cisreg[snp]) for snp in cisreg))

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
	rank = dict((score, index) for index, score in enumerate(sorted(gene_scores.values(), reverse=True)))

	return [ GeneSNP_Association(
		gene = gene,
		snp = snp,
		cisregulatory_evidence = cisreg[gene],
		regulatory_evidence = reg,
		intermediary_scores = intermediary_scores[gene],
		score = gene_scores[gene],
		rank = rank[gene_scores[gene]] + 1)
	for gene in cisreg ]

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
		for evidence in cisreg[gene] + reg:
			if float(evidence.score) > intermediary_scores[gene][evidence.source]:
				intermediary_scores[gene][evidence.source] = float(evidence.score)

			# VEP stats
			if evidence.source == 'VEP':
				intermediary_scores[gene]['VEP_count'] += 1
				intermediary_scores[gene]['VEP_sum'] += float(evidence.score) 

		# Ad hoc bounds defined here:
		# PCHiC
		intermediary_scores[gene]['PCHiC'] = min(intermediary_scores[gene]['PCHiC'], 1)

		# VEP
		if 'VEP' in intermediary_scores[gene]:
			intermediary_scores[gene]['VEP_mean'] = intermediary_scores[gene]['VEP_sum'] / intermediary_scores[gene]['VEP_count']

		# Weighted sum
		gene_scores[gene] = sum(intermediary_scores[gene][source] * postgap.Globals.EVIDENCE_WEIGHTS[source] for source in intermediary_scores[gene] if source in postgap.Globals.EVIDENCE_WEIGHTS)

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
		logging.info("Searching for cis-regulatory data on %i SNPs in all databases" % (len(ld_snps)))
		evidence = concatenate(source().run(ld_snps, tissues) for source in postgap.Cisreg.sources)
	else:
		logging.info("Searching for cis-regulatory data on %i SNPs in (%s)" % (len(ld_snps), ", ".join(postgap.Globals.Cisreg_adaptors)))
		evidence = concatenate(source().run(ld_snps, tissues) for source in postgap.Cisreg.get_filtered_subclasses(postgap.Globals.Cisreg_adaptors))

	filtered_evidence = filter(lambda association: association.gene is not None and association.gene.biotype == "protein_coding", evidence_list)

	# Group by snp, then gene:
	res = collections.defaultdict(lambda: collections.defaultdict(list))
	for association in filtered_evidence:
		assert type(association) is postgap.DataModel.Cisregulatory_Evidence, "association is Cisregulatory_Evidence"
		res[association.snp][association.gene].append(association)

	if postgap.Globals.DEBUG:
		if postgap.Globals.Cisreg_adaptors == None:
			logging.info(("Found %i cis-regulatory interactions in all databases" % (len(res))))
		else:
			logging.info(("Found %i cis-regulatory interactions in (%s)" % (len(res), ", ".join(postgap.Globals.Cisreg_adaptors))))
	return res

def regulatory_evidence(snps, tissues):
	"""

		Extract regulatory evidence linked to SNPs and stores them in a hash
		* [ SNP ]
		* [ string ]
		Returntype: Hash: SNP => [ Regulatory_evidence ]

	"""
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

def gene_to_phenotypes(gene):
	"""

		Look up phenotype annotations for gene
		Args:
		* ENSG stable ID
		Return type: [ OMIM Phenotype ]

	"""
	return [] # TODO and remove stopper
	if gene not in phenotype_cache:
		phenotype_cache[gene['stable_id']] = postgap.Ensembl_lookup.get_gene_phenotypes(gene)

	return phenotype_cache[gene['stable_id']]
