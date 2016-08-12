import sys
import re
import requests
import cPickle as pickle
import collections
import pdb

BIN_WIDTH = 10000;
MAX_DISTANCE = 500000
MAX_BIN = int((MAX_DISTANCE - 1) / BIN_WIDTH)
TSS_cache = dict()
SPECIES = 'Human'
FDR_Model = collections.namedtuple('FDR_Model', ['FDR','BIN_WIDTH','MAX_DISTANCE'])

def main():
	distribution = reduce(process_line, sys.stdin, collections.defaultdict(lambda:  0))
	res = FDR_Model(
		FDR = dict((bin_num, distribution[MAX_BIN] * .99 / distribution[bin_num]) for bin_num in distribution.keys()),
		MAX_DISTANCE = MAX_DISTANCE,
		BIN_WIDTH = BIN_WIDTH
	)
	pickle.dump(res, sys.stdout)

def process_line(distribution, line):
    columns = line.split()[0].split('t') + line.split()[1:]
    chr = columns[0]

    if chr[0] == '#':
        return distribution

    chr = re.sub('^chr','', chr)
    start = int(columns[1])
    end = int(columns[2])
    gene = columns[3]

    if gene not in TSS_cache:
        TSS_cache[gene] = get_gene_tss(gene)
    coords = TSS_cache[gene]
    if coords is None:
        return distribution
    pos = (start + end) / 2
    if coords[0] != chr:
        return distribution
    distance = abs(pos - coords[1])
    if distance > MAX_DISTANCE:
        return distribution
    distribution[int(distance / BIN_WIDTH)] += 1
    return distribution

def get_gene_tss(gene):
	server = "http://grch37.rest.ensembl.org"
	ext = "/lookup/symbol/%s/%s?expand=1" % (SPECIES, gene)

	r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

	if not r.ok:
		return None

	gene = r.json()
	if gene['strand']:
		return (gene['seq_region_name'], gene['start'])
	else:
		return (gene['seq_region_name'], gene['end'])

main()
