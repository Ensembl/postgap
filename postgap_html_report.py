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


def main():

	html_template = 'geneReport.html'

	#get options
	options = get_options()

	#assert that files exists
	assert os.path.exists(options.result_file), "result file " + options.result_file + " can't be found."
	assert os.path.exists(html_template), "HTML template file " + html_template + " can't be found."
	
	#get top 10 genes
	top_10_genes, top_10_snps = get_top_10_genes_snps(options.result_file)

	#get top 10 pathways
	top_10_pathways = get_top_10_pathways(top_10_genes)

	#load template and render the html file out
	template_file = open(html_template, 'r')
	template_data = str(template_file.read())
	t = Template(template_data)
	output_html=open(options.html, "w+") 
	output_html.write(t.render(gene_list=top_10_genes, snp_list=top_10_snps, pathway_list=top_10_pathways))
	output_html.close

def get_options():

	parser = argparse.ArgumentParser(description="Creates a HTML report based on POSTGAP (--output2) results", formatter_class = RawTextHelpFormatter)

	requiredNamed = parser.add_argument_group('required arguments')
	requiredNamed.add_argument('--html', help='Name of the HTML output file', required=True)
	requiredNamed.add_argument('--result_file', help='Name of the POSTGAP (--output2) results file', required=True)
	options = parser.parse_args()

	return options

def get_top_10_genes_snps(result_file):

	top_10_genes = [[0 for x in range(4)] for y in range(10)]
	top_10_snps = [[0 for x in range(4)] for y in range(10)]

	source_file = open(result_file, 'r')
	reader = csv.reader(source_file, delimiter='\t')

	for row in reader:

		gene = row[0]
		cluster = row[1].replace('GWAS_Cluster_', '')
		snp = row[2]
		snp_gene_tissue_posterior = row[3]
		tissue = row[4]
		gene_cluster_tissue_posterior = row[5]

		for i in range(10):
			if top_10_genes[i][0] != gene or top_10_genes[i][1] != cluster or top_10_genes[i][2] != tissue:
				if top_10_genes[i][3] < gene_cluster_tissue_posterior:
					item = [gene, cluster, tissue, gene_cluster_tissue_posterior]
					top_10_genes.insert(i, item)
					del top_10_genes[-1]
					break

		for i in range(10):
			if top_10_snps[i][0] != snp or top_10_snps[i][1] != gene or top_10_snps[i][2] != tissue:
				if top_10_snps[i][3] < snp_gene_tissue_posterior:
					item = [snp, gene, tissue, snp_gene_tissue_posterior]
					top_10_snps.insert(i, item)
					del top_10_snps[-1]
					break

	return top_10_genes, top_10_snps

def get_top_10_pathways(top_10_genes):

	#get the list of genes
	gene_list = []
	for gene in top_10_genes:
		if not gene[0] in gene_list:
			gene_list.append(gene[0])

	#call reactome
	genes_to_analyse = '\n'.join([str(i) for i in gene_list])
	headers = {'Content-Type': 'text/plain',}
	res = requests.post('https://reactome.org/AnalysisService/identifiers/projection/?pageSize=10&page=1', headers=headers, data=genes_to_analyse)
	pathways_jdata = res.json()

	top_10_pathways = []
	for i in range(10):
		pathway_stdi = pathways_jdata['pathways'][i]['stId']
		pathway_name = pathways_jdata['pathways'][i]['name']
		pathway_score = 1 - float(pathways_jdata['pathways'][i]['entities']['fdr'])
		item = [pathway_stdi, pathway_name, pathway_score]
		top_10_pathways.insert(i, item)

	return top_10_pathways


if __name__ == "__main__":
	main()