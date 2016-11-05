import sys
import re
import requests
import cPickle as pickle
import collections
import json

BIN_WIDTH = 10000;
MAX_DISTANCE = 500000
MAX_BIN = int((MAX_DISTANCE - 1) / BIN_WIDTH)
TSS_cache = None
SPECIES = 'Human'

def main():
    lines = sys.stdin.readlines()
    global TSS_cache
    TSS_cache = get_TSSs(lines)
    if TSS_cache is None:
	    sys.exit()
    distribution = reduce(process_line, lines, collections.defaultdict(lambda:  0))

    res = dict(
            [('FDR', dict((bin_num, distribution[MAX_BIN] * .99 / distribution[bin_num]) for bin_num in distribution.keys())),
            ('MAX_DISTANCE', MAX_DISTANCE),
            ('BIN_WIDTH', BIN_WIDTH)]

    )
    pickle.dump(res, sys.stdout)

def get_TSSs(lines):
    genes = dict((gene, 1) for gene in get_gene_symbols(lines))
    return get_all_gene_TSSs(genes)

def get_gene_symbols(lines):
    return map(get_gene_symbol, lines)

def get_gene_symbol(line):
    return line.split()[3]

def process_line(distribution, line):
    columns = line.split('\t')
    chrom = columns[0]

    if chrom[0] == '#':
        return distribution

    start = int(columns[1])
    end = int(columns[2])
    gene = columns[3]

    # Gene was not identified in the first pass
    if gene not in TSS_cache:
	return distribution

    coords = TSS_cache[gene]
    if coords is None:
        return distribution
    pos = (start + end) / 2
    if coords[0] != chrom:
        return distribution
    distance = abs(pos - coords[1])
    if distance >= MAX_DISTANCE:
        return distribution
    distribution[int(distance / BIN_WIDTH)] += 1
    
    return distribution

def get_all_gene_TSSs(genes):
	return dict(concatenate(get_gene_chunk_TSSs(chunk) for chunk in chunks(genes.keys(), 1000)))

def chunks(l, n):
	for i in range(0, len(l), n):
		yield l[i:i+n]

def concatenate(list):
	"""

		Shorthand to concatenate a list of lists
		Args: [[]]
		Returntype: []

	"""
	return sum(filter(lambda elem: elem is not None, list), [])

def get_gene_chunk_TSSs(genes):
	server = "http://grch37.rest.ensembl.org"
	ext = "/lookup/symbol/%s" % (SPECIES)
	data = dict()
	data['symbols'] = genes
	result_hash = get_rest_json(server, ext, data)
	return [get_gene_TSS(gene, result_hash) for gene in result_hash]

def get_gene_TSS(gene, result_hash):
	gene_hash = result_hash[gene]
	if gene_hash['strand']:
		return (gene, (gene_hash['seq_region_name'], gene_hash['start']))
	else:
		return (gene, (gene_hash['seq_region_name'], gene_hash['end']))

def get_rest_json(server, ext, data=None):
	"""
		Args:
		* String (server name)
		* String (extension string)
		Return type: JSON object

	"""
	retries = 0

	while True:
		if data is None:
			headers = { "Content-Type" : "application/json" }
			r = requests.get(str(server)+str(ext), headers = headers)
		else:
			headers = {'Content-Type': 'application/json', 'Accept': 'application/json'}
			r = requests.post(str(server)+str(ext), headers = headers, data = json.dumps(data))

		sys.stderr.write("REST JSON Query: %s%s\n" % (server, ext))

		if not r.ok:
			sys.stderr.write("Failed to get proper response to query %s%s\n" % (server, ext) )
			sys.stderr.write("With headers:\n" + repr(headers) + "\n")
			if data is not None:
				sys.stderr.write("With data:\n" + repr(data) + "\n")
			if 'Retry-After' in r.headers:
				time.sleep(int(r.headers['Retry-After']))
				retries += 1
				continue
			r.raise_for_status()
			sys.exit()

		try:
			return r.json()
		except:
			sys.stderr.write("Failed to get proper response to query %s%s\n" % (server, ext) )
			raise

	# Failed too many times
	sys.exit()

main()
