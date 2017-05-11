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
import sqlite3
import re
import logging
import logging.config

import postgap
import postgap.GWAS
import postgap.Cisreg
import postgap.Reg
import postgap.EFO
import postgap.Globals
import postgap.Integration
from postgap.Utils import *

import sys
from pprint import pformat

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
	
	logging.config.fileConfig('configuration/logging.conf')
	logger = logging.getLogger(__name__)

	options = get_options()
	
	logger.info("Starting postagp with the following options:")
	logger.info(pformat(options))

	if len(options.diseases) > 0 or len(options.efos) > 0:
		res = postgap.Integration.diseases_to_genes(options.diseases, options.efos, "CEPH", options.tissues)
	elif options.rsID is not None:
		res = postgap.Integration.rsIDs_to_genes(options.rsID, options.tissues)
	elif options.coords is not None:
		snp = postgap.DataModel.SNP(rsID = options.coords[0], chrom = options.coords[1], pos = int(options.coords[2]))
		if options.tissues is None:
			options.tissues = ["Whole_Blood"]
		res = postgap.Integration.ld_snps_to_genes([snp], options.tissues)

	if options.db is None:
		if options.output is None:
			output = sys.stdout
		else:
			output = open(options.output, "w")

		if options.json_output:
			formatted_results = json.dumps(objectToDict(res))
		elif options.rsID is None and options.coords is None:
			formatted_results = pretty_output(res)
		else:
			formatted_results = pretty_snp_output(res)

		output.write(formatted_results + "\n")
	else:
		db_output(options.db, res)

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
    parser.add_argument('--coords', nargs=3)
    # parser.add_argument('--populations', nargs='*', default=['1000GENOMES:phase_3:GBR'])
    parser.add_argument('--tissues', nargs='*')
    parser.add_argument('--output')
    parser.add_argument('--db')
    parser.add_argument('--species', nargs='*', default = 'Human')
    parser.add_argument('--database_dir', dest = 'databases', default = 'databases')
    parser.add_argument('--debug', '-g', action = 'store_true')
    parser.add_argument('--json_output', '-j', action = 'store_true')
    options = parser.parse_args()

    postgap.Globals.DATABASES_DIR = options.databases
    postgap.Globals.SPECIES = options.species
    postgap.Globals.DEBUG = postgap.Globals.DEBUG or options.debug

    assert postgap.Globals.DATABASES_DIR is not None
    assert options.rsID is None or (options.efos is None and options.diseases is None)
    assert options.rsID is not None or options.efos is not None or options.diseases is not None or options.coords is not None

    if options.diseases is None:
        options.diseases = []

    if options.efos is None:
        options.efos = filter(lambda X: X is not None, (postgap.EFO.suggest(disease) for disease in options.diseases))

    # Expand list of EFOs to children, concatenate, remove duplicates
    options.efos = concatenate(map(postgap.EFO.children, options.efos))

    return options

def pretty_snp_output(associations):
	"""

		Prints association stats in roughly the same format as STOPGAP for a cluster of SNPs
		Args: [ GeneSNP_Association ]
		Returntype: String

	"""
	header = "\t".join(['snp_rsID', 'chrom', 'pos', 'gene_symbol', 'gene_id', 'score'] + [obj.display_name for obj in postgap.Cisreg.sources + postgap.Reg.sources])
	content = map(pretty_snp_association, associations)
	return "\n".join([header] + content) + "\n"

def pretty_snp_association(association):
	"""

		Prints association stats in roughly the same format as STOPGAP for a cluster of SNPs
		Args: GeneSNP_Association
		Returntype: String

	"""
	snp = association.snp
	gene_name = association.gene.name
	gene_id = association.gene.id
	score = association.score

	results = [snp.rsID, snp.chrom, str(snp.pos), gene_name, gene_id, str(score)]
	results += [str(association.intermediary_scores[functional_source.display_name]) for functional_source in postgap.Cisreg.sources]
	return "\t".join(results)

def pretty_output(associations):
	"""

		Prints association stats in roughly the same format as STOPGAP for a cluster of SNPs
		Args: [ GeneCluster_Association ]
		Returntype: String

	"""
	header = "\t".join(['ld_snp_rsID', 'chrom', 'pos', 'gene_symbol', 'gene_id', 'gene_chrom', 'gene_tss', 'disease_names', 'disease_efo_ids', 'score', 'rank', 'r2', 'gwas_snp_ids', 'ls_snp_is_gwas_snp', 'vep_terms'] + [source.display_name for source in postgap.GWAS.sources + postgap.Cisreg.sources + postgap.Reg.sources])
	content = map(pretty_cluster_association, associations)
	return "\n".join([header] + content)

def pretty_cluster_association(association):
	"""

		Prints association stats in roughly the same format as STOPGAP for a cluster of SNPs
		Args: GeneCluster_Association
		Returntype: String

	"""
	results = genecluster_association_table(association)
	return "\n".join("\t".join(map(unicode, row)) for row in results)

def db_output(db, associations):
	"""

		Writes association stats into DB 
		Args1: string, DB name
		Args2: [ GeneCluster_Association ]

	"""
	conn = sqlite3.connect(db)
	table_sql = '''CREATE TABLE IF NOT EXISTS results (
		ld_snp_rsID TEXT, 
		chrom TEXT, 
		pos INT,
		gene_symbol TEXT,
		gene_id TEXT,
		gene_chrom TEXT,
		gene_tss INT,
		disease_names TEXT,
		disease_efo_ids TEXT, 
		score REAL,
		gwas_snp_ids TEXT,
		ls_snp_is_gwas_snp INT,
		vep_terms TEXT,''' + ",".join([re.sub(" ", "_", source.display_name) + " INT\n" for source in postgap.GWAS.sources + postgap.Cisreg.sources + postgap.Reg.sources]) + ")"
	conn.execute(table_sql)
	map(lambda association: db_output_association(conn, association), associations)
	conn.commit()
	conn.close()

def db_output_association(conn, association):
	"""

		Prints association stats in roughly the same format as STOPGAP for a cluster of SNPs
		Arg1: SQLITE3 connection 
		Arg2: GeneCluster_Association
		Returntype: String

	"""
	results = genecluster_association_table(association)
	sql_cmd = "INSERT INTO results VALUES (%s)" % ",".join("?" for item in results)
	conn.execute(sql_cmd, results)

def genecluster_association_table(association):
	"""

		Returns association stats in roughly the same format as STOPGAP for a cluster of SNPs
		Arg1: GeneCluster_Association
		Returntype: [[ string or float ]]

	"""
	gene_name = association.gene.name
	gene_id = association.gene.id
	gene_chrom = association.gene.chrom
	gene_tss = association.gene.tss
	cluster = association.cluster
	gwas_snps = cluster.gwas_snps
	gwas_snp_rsIDs = [gwas_snp.snp.rsID for gwas_snp in gwas_snps]
	disease_names = list(set(gwas_association.disease.name for gwas_snp in gwas_snps for gwas_association in gwas_snp.evidence))
	disease_efos = list(set(gwas_association.disease.efo for gwas_snp in gwas_snps for gwas_association in gwas_snp.evidence))

	gwas_scores = collections.defaultdict(lambda: collections.defaultdict(lambda: 1))
	for gwas_snp in gwas_snps:
		for gwas_association in gwas_snp.evidence:
			if gwas_association.pvalue < gwas_scores[gwas_association.source][gwas_snp.snp.rsID]:
				gwas_scores[gwas_association.source][gwas_snp.snp.rsID] = gwas_association.pvalue

	results = []
	for gene_snp_association in association.evidence:
		vep_terms = "N/A"
		for evidence in gene_snp_association.cisregulatory_evidence:
			if evidence.source == "VEP":
				vep_terms = ",".join(evidence.info['consequence_terms'])
				break

		row = [
				gene_snp_association.snp.rsID, 
				gene_snp_association.snp.chrom, 
				gene_snp_association.snp.pos, 
				gene_name, 
				gene_id, 
				gene_chrom, 
				gene_tss, 
				",".join(disease_names), 
				",".join(disease_efos), 
				gene_snp_association.score, 
				gene_snp_association.rank,
				association.r2,
				",".join(gwas_snp.snp.rsID for gwas_snp in gwas_snps),
				int(gene_snp_association.snp.rsID in gwas_snp_rsIDs),
				vep_terms
			]

		for gwas_source in postgap.GWAS.sources:
			row.append(",".join(str(gwas_scores[gwas_source.display_name][gwas_snp.snp.rsID]) for gwas_snp in cluster.gwas_snps))
		row += [gene_snp_association.intermediary_scores[functional_source.display_name] for functional_source in postgap.Cisreg.sources + postgap.Reg.sources]
		results.append(row)

	return results

if __name__ == "__main__":
	main()
