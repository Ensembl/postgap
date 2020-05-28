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
parser.add_argument('--submit_interval', type=int, default=90, help='sleeping time (s) added to start the analysis of the next chromosome')
parser.add_argument('--check_interval', type=int, default=300, help='sleeping time (s) added to check again if all the analyses of chromosomes have been done')

args = parser.parse_args()
database_dir = args.database_dir
sumstats = args.summary_stats
submitgap = args.submit_interval
checkgap = args.check_interval

# load functions
def get_checkline(filepath='dir/to/file.tsv'):
	f = open(filepath)
	checkline = f.readline()
	f.close()

	return checkline

def submit_jobs(submitlist, submitgap=90, memory=str(2)):
	"""
	database_dir, tempdir, fnstart, filetype, starttime, checkline take global varables
	"""
	if not all([os.path.exists(tempdir + fnstart + '-chr' + n + '.' + filetype) for n in submitlist]):
		# get which colnumn is variant_id
		ncol_rsid = str(checkline.split().index('variant_id') + 1)
		print('column ' + ncol_rsid + ' is SNP rsid, start splitting by chromosomes and running POSTGAP:')

	for n in submitlist:
		# add a sleeping time between jobs
		if n is not submitlist[0]: 
			print('... wait for ' + str(submitgap) + 's, then submit the next job ...')
			time.sleep(submitgap)

		chrname = fnstart + '-chr' + n
		chrfile = fnstart + '-chr' + n + '.' + filetype
		
		jobname = fnstart + '-chr' + n + '_' + starttime

		# create a new file for each chromosome, if it hasn't been done
		if not os.path.exists(tempdir + chrfile):
			os.system("(head -n1 " + tempdir + filename + " && awk 'NR==FNR {a[$3]; next} $" + ncol_rsid + " in a {print $0}' yalan/variants/variants_rsID_coords-chr" + n + ".tsv " + tempdir + filename + ") > " + tempdir + chrfile)
		else:
			print(chrfile + ' exists in ' + tempdir)

		# res & op2 may exist when re-run jobs
		if os.path.exists(tempdir + chrname + '_res.txt'): os.remove(tempdir + chrname + '_res.txt')
		if os.path.exists(tempdir + chrname + '_op2.txt'): os.remove(tempdir + chrname + '_op2.txt')
		
		# run postgap.py ||| note: -m, -d, results & output2!!!!!!!!!!!
		os.system("gsub -m " + memory + " -n " + jobname + " -d " + tempdir + " -q production-rh74 'python postgap/POSTGAP.py --database_dir " + database_dir + " --summary_stats " + tempdir + chrfile + " --hdf5 /nfs/production/panda/ensembl/funcgen/eqtl/GTEx.V6.88_38.cis.eqtls.h5 --sqlite /nfs/production/panda/ensembl/funcgen/eqtl/GTEx.V6.88_38.cis.eqtls.h5.sqlite3 --bayesian --output " + tempdir + chrname + "_res.txt --output2 " + tempdir + chrname + "_op2.txt -g' -r")
		
		print(time.strftime("%H:%M:%S", time.localtime()) + ', chr' + n + ' is submitted!')

def check_jobs_status(checklist, checkgap=300):
	"""
	tempdir, fnstart, starttime take global varables
	"""
	jobstalist = [False] * len(checklist)

	while not all(jobstalist):
		# check all unfinished jobs
		for n in [n for n,s in zip(checklist, jobstalist) if not s]:
			jobname = fnstart + '-chr' + n + '_' + starttime
			os.system('bjobs -J ' + jobname + ' | cat > ' + tempdir + 'jobstatus.txt')
			#jobstatus.txt is an empty file, which means job isn't running
			jobstalist[checklist.index(n)] = True if os.path.getsize(tempdir + 'jobstatus.txt') == 0 else False

		if all(jobstalist):
			#a job isn't running could be either finished or killed
			print('no job is running now, check if they all quitted normally')
			os.remove(tempdir + 'jobstatus.txt')
			break
		else:
			if int(starttime[-4:]) - checkgap / 3 * 5 < int(time.strftime("%M%S", time.localtime())) < int(starttime[-4:]) + checkgap / 3 * 5:
				print('/'.join([('chr' + n) for n,s in zip(checklist, jobstalist) if not s]) + ' unfinished, check every ' + str(checkgap) + 's')
			time.sleep(checkgap)

def read_job_sum_file(path='dir/', whichjob='jobname', sumtype='out'):
	f = open(path + whichjob + '.' + sumtype)
	text = f.readlines()
	f.close()

	return text

def check_a_chromosome_success(n=str(1)):
	"""
	tempdir, fnstart, starttime take global varables
	"""
	jobname = fnstart + '-chr' + n + '_' + starttime

	outtext = read_job_sum_file(path=tempdir, whichjob=jobname, sumtype='out')

	if 'Successfully completed' in outtext[-23]: complete = True
	elif 'Exited with exit code' in outtext[-23]: complete = False
	else: complete = float('NaN') # a NaN value

	return complete

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
checkline = get_checkline(filepath=tempdir + filename)

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
			checkline = get_checkline(filepath=tempdir + filename)

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

submit_jobs(submitlist=chrlist, submitgap=submitgap)

#5. wait and check if all the jobs are done
time.sleep(300)
print('check whether all the jobs have been done ...')
check_jobs_status(checklist=chrlist, checkgap=checkgap)

#6. check if any has errors, and track running/CPU time
sumdf = pd.DataFrame(columns=['chr', 'success', 'run_time', 'cpu_time', 'max_mem', 'ave_mem'], index=range(len(chrlist)))
sumdf['chr'] = chrlist
for i,n in enumerate(chrlist): sumdf['success'][i] = check_a_chromosome_success(n)

reruntimes = 0
increase_mem = False
while not all(sumdf['success']):
	if reruntimes == 3:
		print('already re-run for 3 times, some consistent error occurs!')
		sys.exit('consistent error occurs, please debug')
	
	if sumdf['success'].isnull().values.any():
		print('/'.join(('chr' + n) for n in sumdf.loc[sumdf['success'].isnull()]['chr']) + ' has an unknown situation!')
		sys.exit('unknown situation occurs, please debug')
	
	rerunlist = []
	for n in sumdf.loc[sumdf['success'] == False]['chr']:
		jobname = fnstart + '-chr' + n + '_' + starttime

		errtext = read_job_sum_file(path=tempdir, whichjob=jobname, sumtype='err')

		if 'HTTPError' in errtext[-1]:
			print('a connection error happened, will re-run chr' + n)
			rerunlist.append(n)
		elif 'KeyboardInterrupt\n' in errtext:
			outtext = read_job_sum_file(path=tempdir, whichjob=jobname, sumtype='out')

			# get Delta Memory
			if float(outtext[-15].split()[-2]) < 0:
				print('memory limit exceeded, will increase memory and re-run chr' + n)
				rerunlist.append(n)
				increase_mem = True
			else:
				print('warnings, warnings, warnings! chr' + n + ' has an unknown error situation!')
		else:
			print('warnings, warnings, warnings! chr' + n + ' has an unknown error situation!')
	
	# if all failed jobs have unexpected errors, stop and debug
	if len(rerunlist) == 0:
		print('/'.join(('chr' + n) for n in sumdf.loc[sumdf['success'] == False]['chr']) + ' has an unknown error situation! warnings, warnings, warnings!')
		sys.exit('unknown error to debug')
	
	reruntimes += 1
	submit_jobs(submitlist=rerunlist, submitgap=submitgap, memory=str(2 + (reruntimes if increase_mem else 0)))
	print('re-run round ' + str(reruntimes) + ': ' + '/'.join(('chr' + n) for n in rerunlist) + ' starts! allocate ' + str(2 + (reruntimes if increase_mem else 0)) + 'G memory')
	
	time.sleep(300)
	print('check whether all the re-run jobs have been done ...')
	check_jobs_status(checklist=rerunlist, checkgap=checkgap)
	
	# update sumdf['success']
	for i,row in sumdf.loc[sumdf['chr'].isin(rerunlist)].iterrows(): sumdf.at[i, 'success'] = check_a_chromosome_success(row['chr'])

print('congratulations!!! all the jobs are finished successfully! let\'s track relevant parameters and wrap all the results together')

for i,n in enumerate(chrlist):
	jobname = fnstart + '-chr' + n + '_' + starttime
	
	outtext = read_job_sum_file(path=tempdir, whichjob=jobname, sumtype='out')

	sumdf.at[i, 'run_time'] = outtext[-11].split()[-2]
	sumdf.at[i, 'cpu_time'] = outtext[-19].split()[-2]
	sumdf.at[i, 'max_mem'] = outtext[-18].split()[-2] if 'MB' in outtext[-18] else 'N/A'
	sumdf.at[i, 'ave_mem'] = outtext[-17].split()[-2] if 'MB' in outtext[-17] else 'N/A'

sumdf.to_csv(filedir + fnstart + '_jobsum_' + starttime + '.txt', header=True, index=False, sep='\t', na_rep='N/A')

#7. combine all the results into one file
print('all the jobs have been completed successfully, and wrap all the results together')
# results file could have only column names, and pd.concat can combine even empty dfs
resdf = pd.concat([pd.read_csv(f, sep='\t', header=0, low_memory=False) for f in glob.glob(tempdir + fnstart + '-chr*_res.txt')], sort=False, ignore_index=True)
resdf.to_csv(filedir + fnstart + '_results_' + starttime + '.txt', header=True, index=False, sep='\t', na_rep='N/A')

# output2 file could be an empty file, which pd.read_csv cannot read
op2df = pd.concat([pd.read_csv(f, sep='\t', header=None, low_memory=False) for f in glob.glob(tempdir + fnstart + '-chr*_op2.txt') if os.path.getsize(f) > 0], sort=False, ignore_index=True)
op2df.to_csv(filedir + fnstart + '_output2_' + starttime + '.txt', header=['Gene_ID', 'Cluster_description', 'SNP_ID', 'CLPP_at_that_SNP', 'Tissue', 'CLPP_over_the_whole_cluster'], index=False, sep='\t', na_rep='N/A')

# remove all the intermediate files
shutil.rmtree(tempdir)
