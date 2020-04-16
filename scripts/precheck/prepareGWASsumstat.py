# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a script file for preworks before running postgap.py
"""

import argparse
import os
import sys
import time

#note: temporarily run this code under the home directory!!!!!!!!!!
dir_gwas = 'yalan/GWAScatalog/'

# get command line arguments
parser = argparse.ArgumentParser(description='to read command line arguments')
parser.add_argument("--filename", required=True, type=str, help="the filename of the uploaded dataset")

args = parser.parse_args()
filename = args.filename

print('start pre-checking for ' + filename)

# save the first line into to a variable for the following checking step
f = open(dir_gwas + filename)
checkline = f.readline()
f.close()

# get the name and the type of file for renaming new generated GWAS files
fnlist = filename.split('.')
filetype = fnlist[-1] #eg:tsv
fnstart = filename.replace('.' + filetype, '') #eg:dubois_test

#1. check if the delimiter is tab
if len(checkline.split('\t')) > 1:
	print('delimiter checking passed')
else:
	print(filename + ' is not separated by tab')

	seplist=[' ', ',', ';', ':', '|']
	t = 0
	for s in seplist:
		if len(checkline.split(s)) > 1:
			print(filename + ' is separated by ' + ('space' if s == ' ' else s) + ', and change the delimiter to tab!')

			# replace delimiter with tab, and save into a new file
			os.rename(dir_gwas + filename, dir_gwas + 'backup_' + filename)
			os.system('cat ' + dir_gwas + 'backup_' + filename + ' | tr "' + s + '" "\t" > ' + dir_gwas + filename)
			os.remove(dir_gwas + 'backup_' + filename)

			# update the checkline
			f = open(dir_gwas + filename)
			checkline = f.readline()
			f.close()

			break
		else:
			t += 1
	if t == len(seplist):
		print('unknown delimiter, please change the delimiter to tab')
		sys.exit("error with delimiter")

# now the separator is tab for sure
#2. check the very important three colnames: variant_id, beta, p-value
if all([impnames in checkline.split() for impnames in ['variant_id', 'beta', 'p-value']]):
	print('colnames checking passed')
else:
	print('please make sure all the following information is included in your data, and please specify them with assigned column names:\nrsid -> "variant_id"\nbeta -> "beta"\np -> "p-value"')
	sys.exit("error with column names")

#3. split GWAS summary file by chromosomes
chrlist = [str(n) for n in range(1, 23)]
chrlist.append('O')

# get which colnumn is variant_id
ncol_rsid = str(checkline.split().index('variant_id') + 1)
print('column ' + ncol_rsid + ' is SNP rsid, splitting starts:')

# if the filetype of the uploaded GWAS summary dataset is ".out", change it to ".tsv". Otherwise it will conflict with the ".out" file from job summary
if filetype == 'out':
	filetype = 'tsv'
	print('change ' + filename + ' to ' + fnstart + '.' + filetype)

for n in chrlist:
	chrfile = fnstart + "-chr" + n + "." + filetype

	# save the dataset of each chromosome into a new summary file
	os.system("(head -n1 " + dir_gwas + filename + " && awk 'NR==FNR {a[$3]; next} $" + ncol_rsid + " in a {print $0}' yalan/variants/variants_rsID_coords-chr" + n + ".tsv " + dir_gwas + filename + ") > " + dir_gwas + chrfile)
	print('chr' + n + ' is done!')

	# run postgap.py |||||| note: -m, -d, results & output2!!!!!!!!!!!
	jobname = fnstart + "-chr" + n
	if n is not str(1): 
		time.sleep(180)
		print("... wait for 3 mins, then submit the next job ...")
	os.system("gsub -m 2 -n " + jobname + " -d " + dir_gwas + " -q production-rh74 'python postgap/POSTGAP.py --database_dir yalan/databases/ --summary_stats " + dir_gwas + chrfile + " --bayesian --output " + dir_gwas + jobname + "_res.txt --output2 " + dir_gwas + jobname + "_op2.txt -g' -r")



#count = 0
#for line in open(chrfile): count += 1

#count = len(open(filename).readlines())

#gsub -m 10 -n Shrine -d . -q production-rh74 'python yalan/GWAScatalog/test.py --filename Shrine.txt' -r 