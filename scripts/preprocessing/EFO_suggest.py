'''

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
	ext = "/summaries/search?query=%s" % url_term
	r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
	 
	if not r.ok:
		sys.stderr.write(repr(r) + "\n")
		sys.stderr.write("Failure searching for term " + term + " with query " + server + ext + "\n")
		r.raise_for_status()
		sys.exit()

	decoded = r.json()
	 
	'''
		Example search result:
		{
			status: "/api/status/ok",
			result: [
				{
					'mid' => '3804279AF8462F3A01EAEE2589C94781F248F9D7',
					'notable' => {
						'name' => 'disease; EFO_0000311',
						'id' => 'summary_from_disease_to_EFO_0000311'
					},
					'name' => 'cancer',
					'score' => '86.11212'
				}
			]	
		}
	'''
	
	hits = decoded['result']
	filtered_hits = filter(lambda X: re.search('EFO_\d+', X['notable']['name']) is not None, hits)
	if len(filtered_hits) == 0:
		return None
	sorted_hits = sorted(filtered_hits, key = lambda X: X['score']) 
	selected_hit = sorted_hits[-1]
	return re.sub(r'^.*(EFO_[0-9]*)$', r'\1', selected_hit['notable']['name'])

main()
