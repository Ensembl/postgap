import sys

for line in sys.stdin:
	items = line.strip().split('\t')
	efo = items[-1]
	while len(items) < 71:
		items.append("")
	items[-1] = efo
	print "\t".join(items)
