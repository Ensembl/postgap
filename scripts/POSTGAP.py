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
import sys
import argparse
import collections
import json

import GWAS
import Cisreg
import Reg
import EFO
import Globals
import Integration
from Utils import *

'''
Development TODO list:

A. Datasets to integrate:
	-Cisregulatory annotations:
		PCHIC (STOPGAP Scoring: Single cell line: +1, multiple cell lines: +2)

	-Epigenetic activity:
		PhyloP (STOPGAP Scoring: FPR 0-0.6: +2, 0.6-0.85: +1,0.85-1: +0)

B. Code improvements:
	Pathways analysis (Downstream)
	Take into account population composition in LD calcs

C. Model improvements:
	Replace PICS with Bayesian model
	Fine mapping of summary data
	Tissue selection

'''

def main():
	"""

		Reads commandline parameters, prints corresponding associated genes with evidence info

	"""
	options = get_options()
	if options.output is None:
		output = sys.stdout
	else:
		output = open(options.output, "w")

	if options.rsID is None:
		res = Integration.diseases_to_genes(options.diseases, options.efos, "CEPH", options.tissues)
	else:
		res = Integration.rsIDs_to_genes(options.rsID, options.tissues)

	if options.json_output:
		formatted_results = "\n".join(map(json.dumps, res))
	elif options.rsID is None:
		formatted_results = pretty_output(res)
	else:
		formatted_results = pretty_snp_output(res)

	output.write(formatted_results)

def get_options():
    """

        Reads commandline parameters
        Returntype:
            {
                diseases: [ string ],
                populations: [ string ],
                tissues: [ string ],
            }

    """
    parser = argparse.ArgumentParser(description='Search GWAS/Regulatory/Cis-regulatory databases for causal genes')
    parser.add_argument('--efos', nargs='*')
    parser.add_argument('--diseases', nargs='*')
    parser.add_argument('--rsID')
    # parser.add_argument('--populations', nargs='*', default=['1000GENOMES:phase_3:GBR'])
    parser.add_argument('--tissues', nargs='*')
    parser.add_argument('--output')
    parser.add_argument('--species', nargs='*', default = 'Human')
    parser.add_argument('--database_dir', dest = 'databases', default = 'databases')
    parser.add_argument('--debug', '-g', action = 'store_true')
    parser.add_argument('--json_output', '-j', action = 'store_true')
    options = parser.parse_args()

    Globals.DATABASES_DIR = options.databases
    Globals.SPECIES = options.species
    Globals.DEBUG = Globals.DEBUG or options.debug

    assert Globals.DATABASES_DIR is not None
    assert options.rsID is None or (options.efos is None and options.diseases is None)
    assert options.rsID is not None or options.efos is not None or options.diseases is not None

    if options.diseases is None:
        options.diseases = []

    if options.efos is None:
        options.efos = filter(lambda X: X is not None, (EFO.suggest(disease) for disease in options.diseases))

    # Expand list of EFOs to children, concatenate, remove duplicates
    options.efos = concatenate(map(EFO.children, options.efos))

    return options

def pretty_snp_output(associations):
	"""

		Prints association stats in roughly the same format as STOPGAP for a cluster of SNPs
		Args: [ GeneSNP_Association ]
		Returntype: String

	"""
	header = "\t".join(['snp_rsID', 'gene_symbol', 'gene_id', 'score'] + [obj.name for obj in Cisreg.sources + Reg.sources])
	content = map(pretty_snp_association, associations)
	return "\n".join([header] + content) + "\n"

def pretty_snp_association(association):
	"""

		Prints association stats in roughly the same format as STOPGAP for a cluster of SNPs
		Args: GeneCluster_Association
		Returntype: String

	"""
	snp = association.snp
	gene_name = association.gene.name
	gene_id = association.gene.id
	score = association.score

	functional_scores = collections.defaultdict(int)
	for evidence in association.regulatory_evidence + association.cisregulatory_evidence:
		functional_scores[evidence.source] += evidence.score
	
	results = [snp.rsID, gene_name, gene_id, str(score)]
	results += [str(functional_scores[functional_source.display_name]) for functional_source in Cisreg.sources]
	return "\t".join(results)


def pretty_output(associations):
	"""

		Prints association stats in roughly the same format as STOPGAP for a cluster of SNPs
		Args: [ GeneCluster_Association ]
		Returntype: String

	"""
	header = "\t".join(['ld_snp_rsID', 'gene_symbol', 'gene_id', 'disease_names', 'disease_efo_ids', 'score', 'gwas_snp_ids'] + [source.display_name for source in GWAS.sources + Cisreg.sources + Reg.sources])
	content = map(pretty_cluster_association, associations)
	return "\n".join([header] + content)

def pretty_cluster_association(association):
	"""

		Prints association stats in roughly the same format as STOPGAP for a cluster of SNPs
		Args: GeneCluster_Association
		Returntype: String

	"""
	gene_name = association.gene.name
	gene_id = association.gene.id
	cluster = association.cluster
	gwas_snps = cluster.gwas_snps
	disease_names = list(set(gwas_association.disease.name for gwas_snp in gwas_snps for gwas_association in gwas_snp.evidence))
	disease_efos = list(set(gwas_association.disease.efo for gwas_snp in gwas_snps for gwas_association in gwas_snp.evidence))

	gwas_scores = collections.defaultdict(lambda: collections.defaultdict(lambda: 1))
	for gwas_snp in gwas_snps:
		for gwas_association in gwas_snp.evidence:
			if gwas_association.pvalue < gwas_scores[gwas_association.source][gwas_snp.snp.rsID]:
				gwas_scores[gwas_association.source][gwas_snp.snp.rsID] = gwas_association.pvalue

	functional_scores = collections.defaultdict(lambda: collections.defaultdict(int))
	snp_scores = collections.defaultdict(int)
	for gene_snp_association in association.evidence:
		for evidence in gene_snp_association.regulatory_evidence + gene_snp_association.cisregulatory_evidence:
			functional_scores[evidence.snp.rsID][evidence.source] += evidence.score
		snp_scores[gene_snp_association.snp.rsID] = gene_snp_association.score
	
	pretty_strings = []
	for ld_snp in cluster.ld_snps:
		if snp_scores[ld_snp.rsID] > 0:
			results = [ld_snp.rsID, gene_name, gene_id, ",".join(disease_names), ",".join(disease_efos), str(snp_scores[ld_snp.rsID])]
			results.append(",".join(gwas_snp.snp.rsID for gwas_snp in gwas_snps))
			for gwas_source in GWAS.sources:
				results.append(",".join(str(gwas_scores[gwas_source.display_name][gwas_snp.snp.rsID]) for gwas_snp in cluster.gwas_snps))
			results += [str(functional_scores[ld_snp.rsID][functional_source.display_name]) for functional_source in Cisreg.sources + Reg.sources]
			pretty_strings.append("\t".join(results))

	return "\n".join(pretty_strings)

if __name__ == "__main__":
	main()
