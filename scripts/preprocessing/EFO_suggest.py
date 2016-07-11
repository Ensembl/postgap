'''

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

	Please email comments or questions to the public Ensembl
	developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

	Questions may also be sent to the Ensembl help desk at
	<http://www.ensembl.org/Help/Contact>.

'''

import sys
import re
import requests

def main():
	phenotype_column = int(sys.argv[1])
	cache = dict()
	for line in sys.stdin:
		# To comply with Unix standards, the column number is provided as 1-based on the command line
		phenotype = line.strip().split('\t')[phenotype_column - 1]

		if phenotype in cache:
			efo = cache[phenotype]
		else:
			efo = efo_suggest(phenotype)
			cache[phenotype] = efo
			if efo is None:
				sys.stderr.write("Could not find EFO for %s\n" % (phenotype))
			else:
				sys.stderr.write("Phenotype %s mapped to %s\n" % (phenotype, efo))

		if efo is None:
			print "%s\tN/A" % (str(line.strip()))
		else:
			print "%s\t%s" % (str(line.strip()), str(efo))

def efo_suggest(term):
	""" efo_suggest
		
		Find most appropriate EFO term for arbitrary string
		Arg:
		* string
		Returntype: string (EFO ID)

	"""
	server = 'http://www.ebi.ac.uk/spot/zooma/v2/api'
	url_term = re.sub("%", "", term)
	url_term = re.sub(" ", "%20", url_term)
	ext = "/summaries?query=%s" % url_term
	r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
	 
	if not r.ok:
		sys.stderr.write(repr(r) + "\n")
		sys.stderr.write("Failure searching for term " + term + " with query " + server + ext + "\n")
		r.raise_for_status()
		sys.exit()

	try:
		'''
			Example search result:
			[
				{
					"id":"00ECD0E64DA2908961D228BF60BB8960C1778F1F",
					"semanticTags":["http://www.ebi.ac.uk/efo/EFO_0000400"],
					"annotationURIs":["http://rdf.ebi.ac.uk/resource/zooma/efo/ABF4EFC8B5BF4598775DC1F6F3532B5B"],
					"annotationSourceURIs":["http://www.ebi.ac.uk/efo/efo.owl"],
					"annotatedPropertyType":"disease",
					"annotatedPropertyValue":"Diabetes",
					"annotatedPropertyUri":"http://rdf.ebi.ac.uk/resource/zooma/5220510DE1D1CF8EEC186E4735AE9BE0",
					"annotationSummaryTypeName":"disease; EFO_0000400",
					"quality":72.803116,
					"uri":"http://rdf.ebi.ac.uk/resource/zooma/annotation_summary/00ECD0E64DA2908961D228BF60BB8960C1778F1F"
				}
			]
		'''
		
		hits = r.json()
		filtered_hits = filter(lambda X: 
					'annotatedPropertyType' in X 
					and X['annotatedPropertyType'] is not None 
					and re.search('disease', X['annotatedPropertyType']) is not None 
					and 'annotationSummaryTypeName' in X 
					and X['annotationSummaryTypeName'] is not None 
					and re.search('EFO', X['annotationSummaryTypeName']) is not None, 
				hits)
		if len(filtered_hits) == 0:
			return None
		sorted_hits = sorted(filtered_hits, key = lambda X: X['quality']) 
		selected_hit = sorted_hits[-1]
		return re.sub(r'^.*(EFO_[0-9]*).*$', r'\1', selected_hit['annotationSummaryTypeName'])
	except:
		sys.stderr.write(repr(r) + "\n")
		sys.stderr.write("Failure searching for term " + term + " with query " + server + ext + "\n")
		raise

main()
