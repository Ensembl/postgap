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
import collections
import math

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
from postgap.Utils import *
import logging

phenotype_cache = ()
PVALUE_CUTOFF = 1e-4
MAX_GWAS_HITS = 500

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
	logger = logging.getLogger(__name__)
	
	res = filter(lambda X: X.pvalue < PVALUE_CUTOFF, scan_disease_databases(diseases, efos))

	logger.info("Found %i GWAS SNPs associated to diseases (%s) or EFO IDs (%s) after p-value filter (%f)" % (len(res), ", ".join(diseases), ", ".join(efos), PVALUE_CUTOFF))

	# Only return the top hits
	return sorted(res, key = lambda gwas_snp: gwas_snp.pvalue)[:MAX_GWAS_HITS]

def scan_disease_databases(diseases, efos):
	"""

		Associates gwas_snps from a list of diseases
		Args:
		* [ string ] (trait descriptions)
		* [ string ] (trait EFO identifiers)
		Returntype: [ GWAS_SNP ]

	"""
	logger = logging.getLogger(__name__)

	logger.info("Searching for GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in all databases" % (", ".join(diseases), ", ".join(efos)))

	gwas_associations = concatenate(source().run(diseases, efos) for source in postgap.GWAS.sources)
	
	logger.info("Done searching for GWAS SNPs. Found %s gwas associations." % len(gwas_associations) )

	associations_by_snp = dict()
	for gwas_association in gwas_associations:
		# Sanity filter to avoid breaking downstream code
		# Looking at you GWAS DB, "p-value = 0", pshaw!
		if gwas_association.pvalue <= 0:
			continue
		if gwas_association.snp in associations_by_snp:
			record = associations_by_snp[gwas_association.snp]
			record.evidence.append(gwas_association)
			
			# If the association is better (has a smaller pvalue) than the 
			# current lead SNP, then make the current one the lead SNP.
			#
			if record.pvalue > gwas_association.pvalue:
				associations_by_snp[gwas_association.snp] = GWAS_SNP(
					snp = record.snp,
					pvalue = gwas_association.pvalue,
					evidence = record.evidence,
				)
		else:
			associations_by_snp[gwas_association.snp] = GWAS_SNP(
				snp = gwas_association.snp,
				pvalue = gwas_association.pvalue,
				# A GWAS association is its own evidence.
				evidence = [ gwas_association ]
			)

	res = associations_by_snp.values()

	logger.info("Found %i unique GWAS SNPs associated to diseases (%s) or EFO IDs (%s) in all databases" % (len(res), ", ".join(diseases), ", ".join(efos)))

	return res

def gwas_snps_to_genes(gwas_snps, populations, tissue_weights):
	"""

		Associates Genes to gwas_snps of interest
		Args:
		* [ GWAS_Association ]
		* { population_name: scalar (weight) }
		* { tissue_name: scalar (weight) }
		Returntype: [ GeneCluster_Association ]

	"""
	logger = logging.getLogger(__name__)
	
	for gwas_snp in gwas_snps:
		assert type(gwas_snp) is GWAS_SNP
	
	# Must set the tissue settings before separating out the gwas_snps
	if tissue_weights is None:
		tissue_weights = gwas_snps_to_tissue_weights(gwas_snps)

	clusters = cluster_gwas_snps(gwas_snps, populations)
	res = concatenate(cluster_to_genes(cluster, tissue_weights, populations) for cluster in clusters)

	logger.info("\tFound %i genes associated to all clusters" % (len(res)))

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
	for gwas_snp in gwas_snps:
		assert type(gwas_snp) is GWAS_SNP
		
	return ['Whole_Blood'] # See FORGE??

def cluster_gwas_snps(gwas_snps, populations):
	"""

		Bundle together gwas_snps within LD threshold
		* [ GWAS_SNP ]
		* [ Bio::Ensembl::Variation::Population ]
		Returntype: [ GWAS_Cluster ]

	"""
	
	logger = logging.getLogger(__name__)
	
	for gwas_snp in gwas_snps:
		assert type(gwas_snp) is GWAS_SNP
	
	gwas_snp_locations = get_gwas_snp_locations(gwas_snps)

	logger.info("Found %i locations from %i GWAS SNPs" % (len(gwas_snp_locations), len(gwas_snps)))

	# For every gwas snp location, create the preclusters.
	#
	# A precluster is a GWAS_Cluster with the current SNP in the 
	# GWAS_Cluster.gwas_snps list and the SNPs found by linkage disequilibrium
	# in GWAS_Cluster.ld_snps.
	# 
	preclusters = filter (lambda X: X is not None, [ gwas_snp_to_precluster(gwas_snp_location, populations) for gwas_snp_location in gwas_snp_locations ])
	
	# Remove preclusters having a SNP in a blacklisted region
	#
	filtered_preclusters = postgap.RegionFilter.region_filter(preclusters)
	
	# Merge precluster that share one of their GWAS_Cluster.ld_snps.
	#
	clusters = merge_preclusters(filtered_preclusters)
	
	if len(clusters)>0:
		
		from postgap.Globals import finemap_gwas_clusters_directory		
		logger.info("Writing %i clusters from %i GWAS SNP locations to %s" % (len(clusters), len(gwas_snp_locations), finemap_gwas_clusters_directory))
		write_gwas_clusters_to_files(clusters, finemap_gwas_clusters_directory)
	
	logger.info("Found %i clusters from %i GWAS SNP locations" % (len(clusters), len(gwas_snp_locations)))

	return clusters

def write_gwas_clusters_to_files(clusters, directory_name_for_storing_the_gwas_clusters):
		
	import os
	if not os.path.exists(directory_name_for_storing_the_gwas_clusters):
		os.makedirs(directory_name_for_storing_the_gwas_clusters)
	
	list_of_files_created = []
	
	for cluster in clusters:

		cluster_file_name = create_file_name_for_gwas_cluster(cluster)
		
		import os.path
		if os.path.isfile(cluster_file_name):
			raise Exception("File " + cluster_file_name + " already exists!")

		f = open(cluster_file_name, 'w')
		
		import pickle
		pickle.dump(cluster, f)
		f.close
		
		list_of_files_created.append(cluster_file_name)
	
	return list_of_files_created
	
def create_file_name_for_gwas_cluster(gwas_cluster):
	
	production_name = create_production_name_for_gwas_cluster(gwas_cluster)
	
	from postgap.Globals import finemap_gwas_clusters_directory
	cluster_file_name = finemap_gwas_clusters_directory + "/" +  production_name + ".pickle"
	
	return cluster_file_name

def create_directory_name_for_storing_the_eqtl_cluster(gwas_cluster):
	
	assert type(gwas_cluster) is postgap.DataModel.GWAS_Cluster, "The gwas_cluster is a GWAS_Cluster "
	
	production_name = create_production_name_for_gwas_cluster(gwas_cluster)
	
	from postgap.Globals import finemap_gwas_clusters_directory
	directory_name_for_storing_the_eqtl_cluster = finemap_gwas_clusters_directory + "/" +  production_name
	
	return directory_name_for_storing_the_eqtl_cluster
	
def create_production_name_for_gwas_cluster(gwas_cluster):
	"""
		Create a production name for a gwas cluster.
		
		Within a run, no other gwas cluster should be given this name, so 
		this should be usable for generating file or directory names. 
	"""
	assert type(gwas_cluster) is postgap.DataModel.GWAS_Cluster, "The gwas_cluster is a GWAS_Cluster "
	
	gwas_snps = gwas_cluster.gwas_snps
	assert type(gwas_snps) is list, "We have a list."
	assert len(gwas_snps) > 0, "We have a list with something in it."
	
	gwas_snp = gwas_snps[0]
	assert type(gwas_snp) is postgap.DataModel.GWAS_SNP, "The gwas_snp is a GWAS_SNP."
	
	snp = gwas_snp.snp
	assert type(snp) is postgap.DataModel.SNP, "snp is a SNP (and not its dbSNP accession)."
	
	cluster_size = len(gwas_snps) + len(gwas_cluster.ld_snps)
	
	rsID = snp.rsID
	production_name = "gwas_cluster_with_" + str(cluster_size) + "_snps_around_" + rsID
	
	return production_name
	
def gwas_snp_to_precluster(gwas_snp, populations):
    """

        Extract neighbourhood of GWAS snp
        Args:
        * [ GWAS_SNP ]
        * [ string ] (populations)
        Returntype: GWAS_Cluster

    """
    logger = logging.getLogger(__name__)
    
    mapped_ld_snps = postgap.LD.calculate_window(gwas_snp.snp)
    logger.info("Found %i SNPs in the vicinity of %s" % (len(mapped_ld_snps), gwas_snp.snp.rsID))
    return GWAS_Cluster(
		gwas_snps = [ gwas_snp ],
		ld_snps = mapped_ld_snps,
		finemap_posteriors = None
	)

def get_gwas_snp_locations(gwas_snps):
	"""

		Extract locations of GWAS SNP:
		Args
		* GWAS_SNP
		Returntype: [ GWAS_SNP ]

	"""
	original_gwas_snp = dict((gwas_snp.snp, gwas_snp) for gwas_snp in gwas_snps)
	mapped_snps = postgap.Ensembl_lookup.get_snp_locations(original_gwas_snp.keys())
	return [
		GWAS_SNP(
			snp = mapped_snp,
			pvalue = original_gwas_snp[mapped_snp.rsID].pvalue,
			evidence = original_gwas_snp[mapped_snp.rsID].evidence,
		)
		for mapped_snp in mapped_snps
		if mapped_snp.rsID in original_gwas_snp
	]

def merge_preclusters(preclusters):
	"""

		Bundle together preclusters that share one LD snp
		* [ Cluster ]
		Returntype: [ Cluster ]

	"""
	logger = logging.getLogger(__name__)
	
	clusters = list(preclusters)
	
	# A dictionary that maps from snp to merged clusters
	snp_owner = dict()
	for cluster in preclusters:
		for ld_snp in cluster.ld_snps:
			# If this SNP has been seen in a different cluster
			if ld_snp in snp_owner and snp_owner[ld_snp] is not cluster:
				
				# Set other_cluster to that different cluster
				other_cluster = snp_owner[ld_snp]

				# Merge data from current cluster into previous cluster
				merged_gwas_snps = other_cluster.gwas_snps + cluster.gwas_snps
				
				#   Create a unique set of ld snps from the current cluster 
				#   and the other cluster.
				merged_ld_snps = dict((ld_snp.rsID, ld_snp) for ld_snp in cluster.ld_snps + other_cluster.ld_snps).values()
				
				#   Create a new Cluster
				merged_cluster = GWAS_Cluster(
					gwas_snps = merged_gwas_snps,
					ld_snps = merged_ld_snps,
					finemap_posteriors = None
				)
				
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

	logger.info("\tFound %i clusters from the GWAS peaks" % (len(clusters)))

	return clusters

def write_geneSNP_Associations_to_files(geneSNP_Associations, directory_name_for_storing_the_eqtl_cluster):
	
	tissue_gene_to_cisregulatory_evidence_list = create_tissue_gene_to_cisreg_dictionary_from_geneSNP_Associations(geneSNP_Associations)
	
	list_of_files_created = []
	
	tissues = tissue_gene_to_cisregulatory_evidence_list.keys()
	for tissue in tissues:
		
		gene_stable_ids = tissue_gene_to_cisregulatory_evidence_list[tissue].keys()		
		for gene_stable_id in gene_stable_ids:
			
			cis_regulatory_evidence_list = tissue_gene_to_cisregulatory_evidence_list[tissue][gene_stable_id]
			assert type(cis_regulatory_evidence_list) is list, "We have a list."
			
			number_of_snps_linked_to_this_gene_and_tissue = len(cis_regulatory_evidence_list)
			
			cluster_file_name = directory_name_for_storing_the_eqtl_cluster + "/cis_regulatory_evidence_from_eqtl_" + str(number_of_snps_linked_to_this_gene_and_tissue) + "_snps_linked_to_" + gene_stable_id + "_in_"  + tissue + ".pickle"
			
			list_of_files_created.append(cluster_file_name)
			
			import os.path
			if os.path.isfile(cluster_file_name):
				raise Exception("File " + cluster_file_name + " already exists!")
			
			import os
			if not os.path.exists(directory_name_for_storing_the_eqtl_cluster):
				os.makedirs(directory_name_for_storing_the_eqtl_cluster)

			f = open(cluster_file_name, 'w')
			import pickle

			for cis_regulatory_evidence in cis_regulatory_evidence_list:
				
				import postgap.DataModel
				assert type(cis_regulatory_evidence) is postgap.DataModel.Cisregulatory_Evidence, "The cis_regulatory_evidence is Cisregulatory_Evidence"
				
				#snp = cis_regulatory_evidence.snp
				#assert type(snp) is postgap.DataModel.SNP, "snp is a SNP (and not its dbSNP accession)."
				
				pickle.dump(cis_regulatory_evidence, f)
				
			f.close
	
	return list_of_files_created
				
def create_tissue_gene_to_cisreg_dictionary_from_geneSNP_Associations(geneSNP_Associations):
	
	import collections
	def generate_default():
		return collections.defaultdict(list)

	tissue_gene_to_cisregulatory_evidence_list = collections.defaultdict(generate_default)
	
	for geneSNP_Association in geneSNP_Associations:
		
		import postgap.DataModel
		assert type(geneSNP_Association) is postgap.DataModel.GeneSNP_Association, "geneSNP_Association is a postgap.DataModel.GeneSNP_Association"
		
		gene = geneSNP_Association.gene
		assert type(gene) is postgap.DataModel.Gene, "gene is a postgap.DataModel.Gene"
		
		gene_stable_id = gene.id
		
		snp = geneSNP_Association.snp
		assert type(snp) is postgap.DataModel.SNP, "snp is a postgap.DataModel.SNP"
		
		cisregulatory_evidence_list = geneSNP_Association.cisregulatory_evidence
		assert type(cisregulatory_evidence_list) is list, "cisregulatory_evidence_list is a list"
		
		for cisregulatory_evidence in cisregulatory_evidence_list:
			
			assert type(cisregulatory_evidence) is postgap.DataModel.Cisregulatory_Evidence, "cisregulatory_evidence is postgap.DataModel.Cisregulatory_Evidence"
			
			tissue = cisregulatory_evidence.tissue
			
			tissue_gene_to_cisregulatory_evidence_list[tissue][gene_stable_id].append(cisregulatory_evidence)
	
	return tissue_gene_to_cisregulatory_evidence_list

def cluster_to_genes(cluster, tissues, populations):
    """

        Associated Genes to a cluster of gwas_snps
        Args:
        * [ Cluster ]
        * { tissue_name: scalar (weights) }
        * { population_name: scalar (weight) }
        Returntype: [ GeneCluster_Association ]

    """
    # Compute LD from top SNP
    top_gwas_hit = sorted(cluster.gwas_snps, key=lambda X: X.pvalue)[-1]
    ld = postgap.LD.get_lds_from_top_gwas(top_gwas_hit.snp, cluster.ld_snps, populations)

    # Obtain interaction data from LD snps
    associations = ld_snps_to_genes([snp for snp in cluster.ld_snps if snp in ld], tissues)
    
    if len(associations) > 0:
	    list_of_files_created = write_geneSNP_Associations_to_files(
			associations, 
			directory_name_for_storing_the_eqtl_cluster = create_directory_name_for_storing_the_eqtl_cluster(cluster)
		)
    
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

    res = [
        GeneCluster_Association(
            gene = gene,
            score = total_score(pics[snp], gene_scores[(gene, snp)][1]),
            cluster = cluster,
            evidence = gene_scores[(gene, snp)][:1], # This is a [ GeneSNP_Association ]
	    r2 = ld[snp]
        )
        for (gene, snp) in gene_scores if snp in pics
    ]

    logger = logging.getLogger(__name__)
    logger.info("\tFound %i genes associated around GWAS SNP %s" % (len(res), top_gwas_hit.snp.rsID))

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
	cisreg = cisregulatory_evidence(ld_snps, tissues)
	
	# Extract SNP specific info:
	reg = regulatory_evidence(cisreg.keys(), tissues)
		
	# Create postgap.DataModel.GeneSNP_Association objects:
	#
	#   cisreg is a dictionary of dictionaries that map to a list of cis regulatory evidence:
	#
	#    dict(postgap.DataModel.SNP -> dict( postgap.DataModel.Gene -> [ postgap.DataModel.Cisregulatory_Evidence ] )
	#
	#   so cisreg[snp] is a dictionary:
	#
	#     dict(postgap.DataModel.Gene -> [ postgap.DataModel.Cisregulatory_Evidence ])
	#
	SNP_GeneSNP_Associations = concatenate((create_SNP_GeneSNP_Associations(snp, reg[snp], cisreg[snp]) for snp in cisreg))
	
	return SNP_GeneSNP_Associations

def find_gene_snp_association_with_source(GeneSNP_Association_list, source):
	
	from postgap.DataModel import GeneSNP_Association
	
	gene_snp_association_with_source = []
	
	for gene_snp_association in GeneSNP_Association_list:
		
		assert type(gene_snp_association) is GeneSNP_Association, "Type is GeneSNP_Association"
		
		if gene_snp_association_has_cisregulatory_evidence_from_source(gene_snp_association, source):
			gene_snp_association_with_source.append(gene_snp_association)

	return gene_snp_association_with_source

def gene_snp_association_has_cisregulatory_evidence_from_source(gene_snp_association, source):
	
	cisregulatory_evidence_list = gene_snp_association.cisregulatory_evidence
	assert type(cisregulatory_evidence_list) is list, "Type is list"
	return cisregulatory_evidence_has_source(cisregulatory_evidence_list, source)
	
	
def cisregulatory_evidence_has_source(cisregulatory_evidence_list, source):
	
	for cisregulatory_evidence in cisregulatory_evidence_list:
		
		assert type(cisregulatory_evidence) is list, "Type is cisregulatory_evidence"
		
		if cisregulatory_evidence.source == source:
			return True
	return False

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
	rank = dict((score, index) for index, score in enumerate(sorted(gene_scores.values())))

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
			intermediary_scores[gene][evidence.source] += float(evidence.score)

			# VEP stats
			if evidence.source == 'VEP':
				if float(evidence.score) > intermediary_scores[gene]['VEP_max']:
					intermediary_scores[gene]['VEP_max'] = float(evidence.score)
				intermediary_scores[gene]['VEP_count'] += 1

		# Ad hoc bounds defined here:
		# PCHiC
		intermediary_scores[gene]['PCHiC'] = min(intermediary_scores[gene]['PCHiC'], 1)

		# VEP
		if 'VEP' in intermediary_scores[gene]:
			intermediary_scores[gene]['VEP_sum'] = intermediary_scores[gene]['VEP']
			intermediary_scores[gene]['VEP'] = intermediary_scores[gene]['VEP_max']
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
	logger = logging.getLogger(__name__)
	
	logger.info("Searching for cis-regulatory data on %i SNPs in all databases" % (len(ld_snps)))
	#evidence = concatenate(source().run(ld_snps, tissues) for source in postgap.Cisreg.sources)
	
	GTEx = postgap.Cisreg.GTEx()
	evidence_list = GTEx.run(ld_snps, tissues)
	
	for evidence in evidence_list:
		assert type(evidence) is postgap.DataModel.Cisregulatory_Evidence, "evidence is Cisregulatory_Evidence"

	filtered_evidence = filter(lambda association: association.gene is not None and association.gene.biotype == "protein_coding", evidence_list)

	# Group by snp, then gene:
	#res = collections.defaultdict(lambda: collections.defaultdict(list))
	res = collections.defaultdict(generate_default)
	for association in filtered_evidence:
		assert type(association) is postgap.DataModel.Cisregulatory_Evidence, "association is Cisregulatory_Evidence"
		res[association.snp][association.gene].append(association)

	logger.info(("Found %i cis-regulatory interactions in all databases" % (len(res))))

	return res

def generate_default():
	return collections.defaultdict(list)

def regulatory_evidence(snps, tissues):
	"""

		Extract regulatory evidence linked to SNPs and stores them in a hash
		* [ SNP ]
		* [ string ]
		Returntype: [ Regulatory_Evidence ]

	"""
	logger = logging.getLogger(__name__)
	
	logger.info("Searching for regulatory data on %i SNPs in all databases" % (len(snps)))

	res = concatenate(source().run(snps, tissues) for source in postgap.Reg.sources)

	logging.info("Found %i regulatory SNPs among %i in all databases" % (len(res), len(snps)))

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

