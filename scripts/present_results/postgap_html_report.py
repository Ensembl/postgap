#! /usr/bin/env python

"""

Copyright [1999-2019] EMBL-European Bioinformatics Institute

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
from argparse import RawTextHelpFormatter
import os
import csv
from jinja2 import Template
import requests
import collections


def main():

	#get options
	options = get_options()

	#assert that files exists
	assert os.path.exists(options.result_file), "result file " + options.result_file + " can't be found."
	assert os.path.exists(options.template), "HTML template file " + options.template + " can't be found."
	
	top_10_genes, top_10_snps, top_10_pathways = get_top_10s(options.result_file)

	#load template and render the html file out
	template_file = open(options.template, 'r')
	template_data = str(template_file.read())
	t = Template(template_data)
	output_html=open(options.output, "w+") 
	output_html.write(t.render(gene_list=top_10_genes, snp_list=top_10_snps, pathway_list=top_10_pathways))
	output_html.close

def get_options():

	parser = argparse.ArgumentParser(description="Creates a HTML report based on POSTGAP (--output2) results", formatter_class = RawTextHelpFormatter)

	requiredNamed = parser.add_argument_group('required arguments')
	requiredNamed.add_argument('--output', help='Name of the HTML output file', required=True)
	requiredNamed.add_argument('--template', help='Name of the HTML template file', required=True)
	requiredNamed.add_argument('--result_file', help='Name of the POSTGAP (--output2) results file', required=True)
	options = parser.parse_args()

	return options

def get_top_10s(result_file):

	gene_cluster_tissue_posterior = collections.defaultdict(float)
	gene_posterior = collections.defaultdict(float)
	snp_table = []

	source_file = open(result_file, 'r')
	reader = csv.reader(source_file, delimiter='\t')

	for row in reader:
		gene, cluster, snp, snp_gene_tissue_posterior, tissue, gene_tissue_posterior = row
		gene_cluster_tissue_posterior[(gene, cluster.replace('GWAS_Cluster_', ''), tissue)] = float(gene_tissue_posterior)
		gene_posterior[gene] += float(gene_tissue_posterior)
		snp_table.append((snp, gene, tissue, snp_gene_tissue_posterior))

	gene_table = [tuple(list(gene_cluster_tissue_tuple) + [gene_cluster_tissue_posterior[gene_cluster_tissue_tuple]]) for gene_cluster_tissue_tuple in gene_cluster_tissue_posterior]
	sorted_gene_table = sorted(gene_table, key=lambda row: -float(row[3]))
	sorted_snp_table = sorted(snp_table, key=lambda row: -float(row[3]))

	#call reactome
	genes_to_analyse = '\n'.join('\t'.join([gene, str(gene_posterior[gene])]) for gene in gene_posterior)
	headers = {'Content-Type': 'text/plain',}
	res = requests.post('https://reactome.org/AnalysisService/identifiers/projection?sortBy=fdr&pageSize=10&page=1', headers=headers, data=genes_to_analyse)
	pathways_jdata = res.json()

	pathway_table = []
	for pathway in pathways_jdata['pathways']:
		pathway_stdi = pathway['stId']
		pathway_name = pathway['name']
		pathway_score = float(pathway['entities']['fdr'])
		pathway_table.append([pathway_stdi, pathway_name, pathway_score])

	sorted_pathway_table = sorted(pathway_table, key=lambda row: row[2])
	return sorted_gene_table[:10], sorted_snp_table[:10], sorted_pathway_table[:10]


if __name__ == "__main__":
	main()
