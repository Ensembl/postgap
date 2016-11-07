import sys

count = int(sys.argv[1])
for line in sys.stdin:
	items = line.strip().split('\t')
	while len(items) < count:
		items.append("")
	print "\t".join(items)
