# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a script file to run postgap.py in parallel
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
parser.add_argument('--memory', type=int, default=2, help='request a node with how large memory (GB) to run a job')

args = parser.parse_args()
database_dir = args.database_dir
sumstats = args.summary_stats
submitgap = args.submit_interval
checkgap = args.check_interval
memory = args.memory

rerunlimts = 3

# load functions
def get_checkline(filepath='dir/to/file.tsv'):
	f = open(filepath)
	checkline = f.readline()
	f.close()

	return checkline

def submit_a_chromosome(chr, memory=str(2)):
	"""
	Note: type(memory) is 'str'
	database_dir, tempdir, fnstart, filetype, starttime take global varables
	"""
	chrname = fnstart + '-chr' + chr
	chrfile = fnstart + '-chr' + chr + '.' + filetype
	jobname = fnstart + '-chr' + chr + '_' + starttime

	# res & op2 may exist when re-run jobs
	if os.path.exists(tempdir + chrname + '_res.txt'): os.remove(tempdir + chrname + '_res.txt')
	if os.path.exists(tempdir + chrname + '_op2.txt'): os.remove(tempdir + chrname + '_op2.txt')
	
	memory = str(memory)
	os.system("gsub -m " + memory + " -n " + jobname + " -d " + tempdir + " -q production-rh74 'python postgap/POSTGAP.py --database_dir " + database_dir + " --summary_stats " + tempdir + chrfile + " --hdf5 /nfs/production/panda/ensembl/funcgen/eqtl/GTEx.V6.88_38.cis.eqtls.h5 --sqlite /nfs/production/panda/ensembl/funcgen/eqtl/GTEx.V6.88_38.cis.eqtls.h5.sqlite3 --bayesian --output " + tempdir + chrname + "_res.txt --output2 " + tempdir + chrname + "_op2.txt -g' -r")
	
	print(time.strftime("%H:%M:%S", time.localtime()) + ', chr' + chr + ' is submitted with ' + memory + 'G memory assigned!')

def start_a_GWAS_sum_stats(chrlist, submitgap=90, memory=str(2)):
	"""
	Note: type(memory) is 'str'
	tempdir, filename, fnstart, filetype, starttime, checkline take global varables
	"""
	# get which colnumn is variant_id
	ncol_rsid = str(checkline.split().index('variant_id') + 1)
	print('column ' + ncol_rsid + ' is SNP rsid, start splitting by chromosomes and running POSTGAP:')

	for n in chrlist:
		# add a sleeping time between jobs
		if n is not chrlist[0]: 
			print('... wait for ' + str(submitgap) + 's, then submit the next job ...')
			time.sleep(submitgap)

		chrname = fnstart + '-chr' + n
		chrfile = fnstart + '-chr' + n + '.' + filetype
		jobname = fnstart + '-chr' + n + '_' + starttime

		# create a new file for each chromosome, and run postgap.py
		os.system("(head -n1 " + tempdir + filename + " && awk 'NR==FNR {a[$3]; next} $" + ncol_rsid + " in a {print $0}' yalan/variants/variants_rsID_coords-chr" + n + ".tsv " + tempdir + filename + ") > " + tempdir + chrfile)
		submit_a_chromosome(chr=n, memory=str(memory))
	
	# remove the original GWAS summary statistics file to save space
	os.remove(tempdir + filename)
	print(filename + ' saved in ' + tempdir + ' has been removed')

def read_job_sum_file(path='dir/', whichjob='jobname', sumtype='out'):
	f = open(path + whichjob + '.' + sumtype)
	text = f.readlines()
	f.close()

	return text

def check_a_chromosome_success(chr):
	"""
	tempdir, fnstart, starttime take global varables
	"""
	jobname = fnstart + '-chr' + chr + '_' + starttime

	success = False
	repeat = False

	# if jobstatus.txt is an empty file, this job isn't running; else, running
	os.system('bjobs -J ' + jobname + ' | cat > ' + tempdir + 'jobstatus.txt')

	# if this job is not running, check whether it finished successfully
	if os.path.getsize(tempdir + 'jobstatus.txt') == 0:
		# just wait for 30s, in case the stdout file has not been completed
		time.sleep(30)

		outtext = read_job_sum_file(path=tempdir, whichjob=jobname, sumtype='out')

		if 'Successfully completed' in outtext[-23]: success = True
		elif 'Exited with exit code' in outtext[-23]: repeat = True
		else: # generate a NaN value, which should never happen
			success = float('NaN')
			print('something should never happen happened to chr ' + chr + '! that line was "' + outtext[-23].rstrip() + '" when checking')
	
	os.remove(tempdir + 'jobstatus.txt')

	return success, repeat

def increse_memory(chr):
	"""
	tempdir, fnstart, starttime take global varables
	"""
	jobname = fnstart + '-chr' + chr + '_' + starttime
	errtext = read_job_sum_file(path=tempdir, whichjob=jobname, sumtype='err')

	increase_mem = False

	if 'HTTPError' in errtext[-1]:
		print('a http error happened, no need to increase memory when re-run chr' + chr)
	elif 'ConnectionError' in errtext[-1]:
		print('a connection error happened, no need to increase memory when re-run chr' + chr)
	elif 'IOError' in errtext[-1]:
		print('chr' + chr + ' gave ' + errtext[-1].rstrip() + '. performance issue with the filesystem, just try again')
	elif 'KeyboardInterrupt\n' in errtext:
		# check stdout file to figure why it was killed
		outtext = read_job_sum_file(path=tempdir, whichjob=jobname, sumtype='out')

		# get Delta Memory
		if float(outtext[-15].split()[-2]) < 0:
			print('memory limit exceeded, will increase memory and re-run chr' + chr)
			increase_mem = True
		else:
			print('warnings, warnings, warnings! chr' + chr + ' has an unknown error situation!')
			sys.exit('unknown error situation, please debug')
	else:
		print('warnings, warnings, warnings! chr' + chr + ' has an unknown error situation!')
		sys.exit('unknown error situation, please debug')

	return increase_mem

def trace_chromosomes_status(chrlist, init_mem=2, checkgap=checkgap, rerunlimts=3):
	"""
	Note: type(init_mem) is 'int'
	tempdir, fnstart, filetype, starttime take global varables
	"""
	sumdf = pd.DataFrame(columns=['chr', 'success', 'repeats', 'run_time', 'cpu_time', 'max_mem', 'ave_mem'], index=range(len(chrlist)))
	sumdf['chr'] = chrlist
	sumdf['repeats'] = [0] * len(chrlist)
	
	while not all(sumdf['success'] == True):
		# check all unfinished jobs
		for i,row in sumdf.loc[sumdf['success'] != True].iterrows():
			chr = row['chr']

			sumdf.at[i, 'success'], repeat = check_a_chromosome_success(chr=chr)

			# when one finished successfully, record relevant parameters
			if sumdf.loc[i, 'success'] is True:
				jobname = fnstart + '-chr' + chr + '_' + starttime
				outtext = read_job_sum_file(path=tempdir, whichjob=jobname, sumtype='out')
				
				sumdf.ix[i, 'run_time':'ave_mem'] = [outtext[-11].split()[-2], outtext[-19].split()[-2], outtext[-18].split()[-2] if 'MB' in outtext[-18] else 'N/A', outtext[-17].split()[-2] if 'MB' in outtext[-17] else 'N/A']

				# once succeeds, remove the chromosomes GWAS summary statistics file to save space
				chrfile = fnstart + '-chr' + chr + '.' + filetype
				os.remove(tempdir + chrfile)

			if repeat:
				sumdf.at[i, 'repeats'] += 1
				reruntimes = sumdf.loc[i, 'repeats']
				
				if reruntimes <= rerunlimts:
					print('chr' + chr + '\'s re-run round ' + str(reruntimes) + ' starts')
					increase_mem = increse_memory(chr=chr)
					submit_a_chromosome(chr=chr, memory=str(round(int(init_mem) * (1 + (reruntimes if increase_mem else 0) / 2))))

		if all(sumdf['success'] == True):
			print('congratulations!!! all the jobs finished successfully!')
			break
		elif any(sumdf['success'].isnull()):
			print('/'.join(('chr' + n) for n in sumdf.loc[sumdf['success'].isnull()]['chr']) + ' quitted with an unknown situation!')
			sys.exit('unknown situation occurs, please debug')
		elif any((sumdf['success'] == False) & (sumdf['repeats'] > rerunlimts)):
			print('/'.join(('chr' + n) for n in sumdf.loc[(sumdf['success'] == False) & (sumdf['repeats'] > rerunlimts), 'chr']) + ' has already re-ran for ' + str(rerunlimts) + ' times, some consistent error occurs!')
			sys.exit('consistent error occurs, please debug')

		if int(starttime[-4:]) - checkgap / 3 * 5 < int(time.strftime("%M%S", time.localtime())) < int(starttime[-4:]) + checkgap / 3 * 5:
			print('/'.join(('chr' + n) for n in sumdf.loc[sumdf['success'] != True, 'chr']) + ' unfinished, check every ' + str(checkgap) + 's')
		time.sleep(checkgap)
	
	return sumdf

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

#4. split GWAS summary file by chromosomes and run postgap.py on each
chrlist = [str(n) for n in range(1, 23)]
chrlist.append('O')

start_a_GWAS_sum_stats(chrlist=chrlist, submitgap=submitgap, memory=str(memory))

#5. keep tracing until all the jobs completed successfully, and make a record of running/CPU time and max/average memory
time.sleep(300)
print('start checking whether all the jobs completed successfully ...')

sumdf = trace_chromosomes_status(chrlist=chrlist, init_mem=memory, checkgap=checkgap, rerunlimts=rerunlimts)
sumdf.to_csv(filedir + fnstart + '_jobsum_' + starttime + '.txt', header=True, index=False, sep='\t', na_rep='N/A')

#6. combine all the results into one file
print('wrap all the results together, and clean up')

# results file could have only column names, and pd.concat can combine even empty dfs
resdf = pd.concat([pd.read_csv(f, sep='\t', header=0, low_memory=False) for f in glob.glob(tempdir + fnstart + '-chr*_res.txt')], sort=False, ignore_index=True)
resdf.to_csv(filedir + fnstart + '_results_' + starttime + '.txt', header=True, index=False, sep='\t', na_rep='N/A')

# output2 file could be an empty file, which pd.read_csv cannot read
op2df = pd.concat([pd.read_csv(f, sep='\t', header=None, low_memory=False) for f in glob.glob(tempdir + fnstart + '-chr*_op2.txt') if os.path.getsize(f) > 0], sort=False, ignore_index=True)
op2df.to_csv(filedir + fnstart + '_output2_' + starttime + '.txt', header=['Gene_ID', 'Cluster_description', 'SNP_ID', 'CLPP_at_that_SNP', 'Tissue', 'CLPP_over_the_whole_cluster'], index=False, sep='\t', na_rep='N/A')

# remove all the intermediate files
for f in glob.glob(tempdir + fnstart + '-chr*_*.txt'): os.remove(f)
#shutil.rmtree(tempdir)
