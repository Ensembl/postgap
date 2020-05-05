# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a script file for preworks before running postgap.py in parallel
"""

import argparse
import glob
import os
import pandas as pd
import shutil
import sys
import time

#note: temporarily run this code under the home directory!!!!!!!!!!
starttime = time.strftime("%y%m%d%H%M%S", time.localtime()) #eg:200501154145

# get command line arguments
parser = argparse.ArgumentParser(description='to read command line arguments')
parser.add_argument('--database_dir', type=str, default='databases/', help='directory where data files are stored')
parser.add_argument('--summary_stats', required=True, type=str, help='location of input summary statistics file')
parser.add_argument('--time_interval', type=int, default=180, help='sleeping time (s) added between analyses of chromosomes')

args = parser.parse_args()
database_dir = args.database_dir
sumstats = args.summary_stats
timegap = args.time_interval

#0. preparations
# get basic info about the GWAS summary statistics
filename = sumstats.split('/')[-1]
filedir = sumstats.replace(filename, '') if len(sumstats.split('/')) > 1 else './'
filetype = filename.split('.')[-1] #eg:tsv
fnstart = filename.replace('.' + filetype, '') #eg:WojcikGL
print(filename + ' is in the ' + ('current directory' if filedir == './' else ('directory ' + filedir)) + '\nstart pre-checking for ' + fnstart)

# create a temporary folder for all the intermediate files that generated during the analysis
tempdir = filedir + fnstart + '_' + starttime + '/'
os.makedirs(tempdir)

#1. check the filetype of the GWAS summary statistics
#a. if it is an excel file, convert it into a tsv file
if filetype in ['xls', 'xlsx', 'xlsm', 'xlsb']:
	pd.read_excel(sumstats).to_csv(tempdir + fnstart + '.tsv', header=True, index=False, sep='\t')
	print(filename + ' is an excel file, convert it into a tab-separated file')
	# update filetype and filename
	filetype = 'tsv'
	filename = fnstart + '.' + filetype
else:
	#b. if it is ".out", change it to ".tsv". Note: conflict with the ".out" file from job summary
	if filetype == 'out':
		os.rename(sumstats, filedir + fnstart + '.tsv')
		print('change filename from ' + filename + ' to ' + fnstart + '.tsv')
		# update filetype and filename
		filetype = 'tsv'
		filename = fnstart + '.' + filetype

	shutil.copyfile(filedir + filename, tempdir + filename)

# save the first line into to a variable for the following checking step
f = open(tempdir + filename)
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
			os.rename(tempdir + filename, tempdir + fnstart + '_backup.' + filetype)
			os.system('cat ' + tempdir + fnstart + '_backup.' + filetype + ' | tr "' + s + '" "\t" > ' + tempdir + filename)

			# update the checkline
			f = open(tempdir + filename)
			checkline = f.readline()
			f.close()

			break
		else:
			t += 1
	if t == len(seplist):
		print('unknown delimiter, please change the delimiter to tab')
		sys.exit('error with delimiter')

# now the separator is tab for sure
#3. check three very important colnames: variant_id, beta, p-value
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
	os.system("(head -n1 " + tempdir + filename + " && awk 'NR==FNR {a[$3]; next} $" + ncol_rsid + " in a {print $0}' yalan/variants/variants_rsID_coords-chr" + n + ".tsv " + tempdir + filename + ") > " + tempdir + chrfile)
	os.system("gsub -m 2 -n " + jobname + '_' + starttime + " -d " + tempdir + " -q production-rh74 'python postgap/POSTGAP.py --database_dir " + database_dir + " --summary_stats " + tempdir + chrfile + " --hdf5 /nfs/production/panda/ensembl/funcgen/eqtl/GTEx.V6.88_38.cis.eqtls.h5 --sqlite /nfs/production/panda/ensembl/funcgen/eqtl/GTEx.V6.88_38.cis.eqtls.h5.sqlite3 --bayesian --output " + tempdir + jobname + "_res.txt --output2 " + tempdir + jobname + "_op2.txt -g' -r")
	
	print(time.strftime("%H:%M:%S", time.localtime()) + ', chr' + n + ' is submitted!')

#5. wait and check if all the jobs are done
#time.sleep(60)
print('check whether all the jobs have been done ...')

jobstalist = [False] * len(chrlist)
while not all(jobstalist):
	# check all unfinished jobs
	for n in [i for i,s in zip(chrlist, jobstalist) if not s]:
		jobname = fnstart + '-chr' + n
		os.system('bjobs -J ' + jobname + '_' + starttime + ' | cat > ' + tempdir + 'jobstatus.txt')
		#jobstatus.txt is an empty file, which means job isn't running
		jobstalist[-1 if n == chrlist[-1] else (int(n) - 1)] = True if os.path.getsize(tempdir + 'jobstatus.txt') == 0 else False

	if all(jobstalist):
		#a job isn't running could be either finished or killed
		print('no job is running now, check if they all quitted normally')
		os.remove(tempdir + 'jobstatus.txt')
		break
	else:
		print('/'.join([('chr' + i) for i,s in zip(chrlist, jobstalist) if not s]) + ' unfinished, wait for 30s and check again ...')
		time.sleep(30)

#6. check if any has errors, and track running/CPU time
sumdf = pd.DataFrame(columns=['chr', 'success', 'run_time', 'cpu_time', 'max_mem', 'ave_mem'], index=range(len(chrlist)))

for i,n in enumerate(chrlist):
	jobname = fnstart + '-chr' + n

	f = open(tempdir + jobname + '_' + starttime + '.out')
	outtext = f.readlines()
	f.close()

	if 'Successfully completed' in outtext[-23]: complete = True
	elif 'Exited with exit code' in outtext[-23]: complete = False
	else: complete = 'N/A'

	sumdf.iloc[i] = [n, complete, outtext[-11].split()[-2], outtext[-19].split()[-2], outtext[-18].split()[-2] if 'MB' in outtext[-18] else 'N/A', outtext[-17].split()[-2] if 'MB' in outtext[-17] else 'N/A']

sumdf.to_csv(filedir + fnstart + '_jobsum_' + starttime + '.txt', header=True, index=False, sep='\t', na_rep='N/A')

#7. combine all the results into one file
if all(sumdf['success']):
	print('all the jobs have been completed successfully, and wrap all the results together')
	# results file could have only column names, and pd.concat can combine even empty dfs
	resdf = pd.concat([pd.read_csv(f, sep='\t', header=0) for f in glob.glob(tempdir + fnstart + '-chr*_res.txt')], sort=False, ignore_index=True)
	resdf.to_csv(filedir + fnstart + '_results_' + starttime + '.txt', header=True, index=False, sep='\t', na_rep='N/A')

	# output2 file could be an empty file, which pd.read_csv cannot read
	op2df = pd.concat([pd.read_csv(f, sep='\t', header=None) for f in glob.glob(tempdir + fnstart + '-chr*_op2.txt') if os.path.getsize(f) > 0], sort=False, ignore_index=True)
	op2df.to_csv(filedir + fnstart + '_output2_' + starttime + '.txt', header=['Gene_ID', 'Cluster_description', 'SNP_ID', 'CLPP_at_that_SNP', 'Tissue', 'CLPP_over_the_whole_cluster'], index=False, sep='\t', na_rep='N/A')

	# remove all the intermediate files
	shutil.rmtree(tempdir)
else:
	print('some gave error warnings, please check!!!')
