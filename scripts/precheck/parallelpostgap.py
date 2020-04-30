# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a script file for preworks before running postgap.py in parallel
"""

import argparse
import os
import sys
import time
import glob
import pandas as pd

#note: temporarily run this code under the home directory!!!!!!!!!!

# get command line arguments
parser = argparse.ArgumentParser(description='to read command line arguments')
parser.add_argument('--database_dir', type=str, default='databases/', help='directory where data files are stored')
parser.add_argument('--summary_stats', required=True, type=str, help='location of input summary statistics file')
parser.add_argument('--time_interval', type=int, default=180, help='sleeping time (s) added between analyses of chromosomes')

args = parser.parse_args()
database_dir = args.database_dir
sumstats = args.summary_stats
timegap = args.time_interval

#1. check the filetype of the GWAS summary statistics
filetype = sumstats.split('.')[-1] #eg:xlsx

#a. if it is an excel file, convert it into a tsv file
if filetype in ['xls', 'xlsx', 'xlsm', 'xlsb']:
	sumstats_old = sumstats
	sumstats = sumstats.replace(filetype, 'tsv')
	pd.read_excel(sumstats_old).to_csv(sumstats, header=True, index=False, sep='\t')
	#os.remove(sumstats_old)
	print(sumstats_old.split('/')[-1] + ' is an excel file, convert it into a tab-separated file')

#b. if it is ".out", change it to ".tsv". Note: conflict with the ".out" file from job summary
if filetype == 'out':
	sumstats_old = sumstats
	sumstats = sumstats.replace(filetype, 'tsv')
	os.rename(sumstats_old, sumstats)
	print('change filename ' + sumstats_old.split('/')[-1] + ' to ' + sumstats.split('/')[-1])

# get relevant info
filename = sumstats.split('/')[-1]
filedir = sumstats.replace(filename, '') if len(sumstats.split('/')) > 1 else './'
print(filename + ' is in the ' + ('current directory' if filedir == './' else ('directory ' + filedir)))
filetype = filename.split('.')[-1] #eg:tsv
fnstart = filename.replace('.' + filetype, '') #eg:dubois_test
print('start pre-checking for ' + fnstart)

# save the first line into to a variable for the following checking step
f = open(sumstats)
checkline = f.readline()
f.close()

#2. check if the delimiter is tab
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
			os.rename(sumstats, filedir + 'backup_' + filename)
			os.system('cat ' + filedir + 'backup_' + filename + ' | tr "' + s + '" "\t" > ' + sumstats)
			os.remove(filedir + 'backup_' + filename)

			# update the checkline
			f = open(filedir + filename)
			checkline = f.readline()
			f.close()

			break
		else:
			t += 1
	if t == len(seplist):
		print('unknown delimiter, please change the delimiter to tab')
		sys.exit('error with delimiter')

# now the separator is tab for sure
#3. check the very important three colnames: variant_id, beta, p-value
if all([impnames in checkline.split() for impnames in ['variant_id', 'beta', 'p-value']]):
	print('colnames checking passed')
else:
	print('please make sure all the following information is included in your data, and please specify them with assigned column names:\nrsid -> "variant_id"\nbeta -> "beta"\np -> "p-value"')
	sys.exit('error with column names')

#4. split GWAS summary file by chromosomes and run postgap.py
chrlist = [str(n) for n in range(1, 23)]
chrlist.append('O')

# get which colnumn is variant_id
ncol_rsid = str(checkline.split().index('variant_id') + 1)
print('column ' + ncol_rsid + ' is SNP rsid, splitting starts:')

for n in chrlist:
	# add a sleeping time between jobs
	if n is not str(1): 
		print('... wait for ' + str(timegap) + 's, then submit the next job ...')
		time.sleep(timegap)

	chrfile = fnstart + '-chr' + n + '.' + filetype
	jobname = fnstart + '-chr' + n

	# create a new file for each chromosome, and run postgap.py |||||| note: -m, -d, results & output2!!!!!!!!!!!
	os.system("(head -n1 " + sumstats + " && awk 'NR==FNR {a[$3]; next} $" + ncol_rsid + " in a {print $0}' yalan/variants/variants_rsID_coords-chr" + n + ".tsv " + sumstats + ") > " + filedir + chrfile)
	os.system("gsub -m 2 -n " + jobname + " -d " + filedir + " -q production-rh74 'python postgap/POSTGAP.py --database_dir " + database_dir + " --summary_stats " + filedir + chrfile + " --hdf5 /nfs/production/panda/ensembl/funcgen/eqtl/GTEx.V6.88_38.cis.eqtls.h5 --sqlite /nfs/production/panda/ensembl/funcgen/eqtl/GTEx.V6.88_38.cis.eqtls.h5.sqlite3 --bayesian --output " + filedir + jobname + "_res.txt --output2 " + filedir + jobname + "_op2.txt -g' -r")
	
	print('chr' + n + ' is running!')

#5. wait and check if all the jobs are done
#time.sleep(120)
print('check whether all the jobs have been done ...')

jobstalist = [False] * len(chrlist)
while not all(jobstalist):
	# check all unfinished jobs
	for n in [i for i,s in zip(chrlist, jobstalist) if not s]:
		jobname = fnstart + '-chr' + n
		os.system('bjobs -J ' + jobname + ' | cat > jobstatus.txt')
		#jobstatus.txt is an empty file, which means job isn't running
		jobstalist[-1 if n == chrlist[-1] else (int(n) - 1)] = True if os.path.getsize('jobstatus.txt') == 0 else False

	if all(jobstalist):
		#a job isn't running could be either finished or killed
		print('no job is running now, check if they all quitted normally')
		os.remove('jobstatus.txt')
		break
	else:
		print(', '.join([('chr' + i) for i,s in zip(chrlist, jobstalist) if not s]) + ' unfinished, wait for 30s and check again ...')
		time.sleep(30)

#6. combine all the results into one file
if all([os.path.getsize(f) == 0 for f in glob.glob(filedir + fnstart + '-chr*.err')]):
	print('no error warnings, wrap all the results together')

	# results file could have only column names, and pd.concat can combine even empty dfs
	resdf = pd.concat([pd.read_csv(f, sep='\t', header=0) for f in glob.glob(filedir + fnstart + '-chr*_res.txt')], sort=False, ignore_index=True)
	resdf.to_csv(filedir + fnstart + '_results.txt', header=True, index=False, sep='\t', na_rep='N/A')

	# output2 file could be an empty file, which pd.read_csv cannot read
	op2df = pd.concat([pd.read_csv(f, sep='\t', header=None) for f in glob.glob(filedir + fnstart + '-chr*_op2.txt') if os.path.getsize(f) > 0], sort=False, ignore_index=True)
	op2df.to_csv(filedir + fnstart + '_output2.txt', header=['Gene_ID', 'Cluster_description', 'SNP_ID', 'CLPP_at_that_SNP', 'Tissue', 'CLPP_over_the_whole_cluster'], index=False, sep='\t', na_rep='N/A')

	# remove all the intermediate files
	for f in glob.glob(filedir + fnstart + '-chr*'): os.remove(f)
else:
	print('some gave error warnings!!!')
