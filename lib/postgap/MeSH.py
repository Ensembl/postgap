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
import sys
import requests
import xmltodict

NCBI_Taxon_ID = {
	'Human': 9606
}

def gene_to_postgap.MeSH(gene):
	"""

		Look up postgap.MeSH annotations for gene
		Args:
		* [ string ] (gene names)
		Return type: [ string ] (annotations)
	"""
	server = "http://gene2mesh.ncibi.org"
	ext = "/fetch?genesymbol=%s&taxid=%s" % (gene.name, NCBI_Taxon_ID[postgap.Globals.SPECIES])
	print '>>>>>>>>>>>>>>>>>>'
	print str(server)+str(ext)
	print '>>>>>>>>>>>>>>>>>>'
	response = requests.get(str(server)+str(ext))

	if not response.ok:
		sys.stderr.write("Failed to get proper response to query %s%s\n" % (server, ext) )
		sys.stderr.write(response.content + "\n")
		response.raise_for_status()
		print repr(response)

	'''
	Example postgap.MeSH output:

	{
		'Gene2postgap.MeSH': {
			'Request': {
				'ParameterSet': {
					'Tool': 'none',
					'GeneSymbol': 'csf1r',
					'TaxonomyID': '9606',
					'Limit': '1000',
					'Email': 'anonymous'
				},
				'type': 'fetch'
			},
			'Response': {
				'ResultSet': {
					'Result': [
						{
							'FisherExact': {
								'content': '1.8531319238671E-230',
								'type': 'p-value'
							},
							'ChiSquare': '112213.6506462',
							'Fover': '1498.1813411401',
							'postgap.MeSH': {
								'Qualifier': {
									'Name': 'metabolism'
								},
								'Descriptor': {
									'TreeNumber': [
										'D08.811.913.696.620.682.725.400.500',
										'D12.776.543.750.060.492',
										'D12.776.543.750.705.852.150.150',
										'D12.776.543.750.750.400.200.200',
										'D12.776.624.664.700.800',
										'D23.050.301.264.035.597',
										'D23.101.100.110.597'
									],
									'Identifier': 'D016186',
									'Name': 'Receptor, Macrophage Colony-Stimulating Factor',
									'UMLSID': {}
								}
							},
							'DocumentSet': {
								'type': 'pubmed',
								'PMID': [
								]
							},
							'Gene': {
								'Taxonomy': {
									'Identifier': '9606'
								},
								'Identifier': '1436',
								'type': 'Entrez',
								'Description': 'colony stimulating factor 1 receptor',
								'Symbol': 'CSF1R'
							}
						},
					],
					'sort': 'FisherExact',
					'count': '94',
					'order': 'ascending'
				},
				'Copyright': {
					'Details': 'http://nlp.ncibi.org/Copyright.txt',
					'Year': '2009',
					'Statement': 'Copyright 2009 by the postgap.Regents of the University of Michigan'
				},
				'Support': {
					'Details': 'http://www.ncibi.org',
					'GrantNumber': 'U54 DA021519',
					'Statement': 'Supported by the National Institutes of Health as part of the NIH\\\'s National Center for Integrative Biomedical Informatics (NCIBI)'
				}
			}
		}
	}
	'''

	print response.content
	if len(response.content):
		try:
			hash = xmltodict.parse(response.content)
			print repr(hash)
			hits = hash['Gene2postgap.MeSH']['Request']['ResultSet']['Result']
			# XML subtletly: if a tag is repeated, a list of objects is produced,
			# else a single object. Careful when unpacking!
			if hits is list:
				return [hit['postgap.MeSH']['Descriptor']['Name'] for hit in hits]
			else:
				return [hits['postgap.MeSH']['Descriptor']['Name']]
		except:
			return []
	else:
		return []

