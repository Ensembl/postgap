import csv
import sys

reader = csv.reader(sys.argv[1])
for row in reader:
	print "\t".join(row)
