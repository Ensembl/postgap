import sys
import cPickle as pickle

file = open('postgap_output')
data = pickle.load(file)
file.close()

for genecluster_association in data:
	gene = genecluster_association.gene
	cluster = genecluster_association.cluster
	chrom = cluster.ld_snps[0].chrom
	start = min(X.pos for X in cluster.ld_snps)
	end = max(X.pos for X in cluster.ld_snps)
	posterior = genecluster_association.collocation_posterior
	for tissue in posterior:
		print "\t".join(map(str, [gene.name, gene.id, tissue, chrom, start, end, posterior[tissue]]))
