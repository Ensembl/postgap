import csv
import sys

reader = csv.reader(sys.stdin)
for row in reader:
	print "\t".join(row)
