'''

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

	Please email comments or questions to the public Ensembl
	developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

	Questions may also be sent to the Ensembl help desk at
	<http://www.ensembl.org/Help/Contact>.

'''

import sys
import re
import requests
import postgap.EFO

phenotype_column = int(sys.argv[1])
cache = dict()

# Preload cache with precomupted suggestions:
for filename in sys.argv[2:]:
	file = open(filename)
	for line in file:
		items = line.strip().split('\t')
		string = items[0]

		while string[0] == '"' and string[-1] == '"':
			string = string[1:-1]
			
		if string in cache:
			cache[string] = cache[string] + "," + items[1]
		else:
			cache[string] = items[1]

for line in sys.stdin:
	# To comply with Unix standards, the column number is provided as 1-based on the command line
	phenotype = str(line.strip().split('\t')[phenotype_column - 1])
	if phenotype[0] == '"' and phenotype[-1] == '"':
		phenotype = phenotype[1:-1]

	if phenotype in cache:
		efo = cache[phenotype]
	else:
		sys.stderr.write(phenotype + "\n")
		try:
			uni = unicode(phenotype)
			efo = postgap.EFO.suggest(uni)
		except UnicodeDecodeError:
			uni = unicode(phenotype.decode('latin1'))		
			efo = postgap.EFO.suggest(uni)

		efo = postgap.EFO.suggest(uni)

		cache[phenotype] = efo
		if efo is None:
			sys.stderr.write("Could not find EFO for %s\n" % (phenotype))
		else:
			sys.stderr.write("Phenotype %s mapped to %s\n" % (phenotype, str(efo)))

	if efo is None:
		print "%s\tN/A" % (str(line[:-1]))
	else:
		print "%s\t%s" % (str(line[:-1]), str(efo))
