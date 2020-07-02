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
parser.add_argument('--memory', type=int, default=2, help='memory of a node (GB) requested to run a job')
parser.add_argument('--kstart', type=int, default=1, help='how many causal variants to start with in the full exploration of sets')
parser.add_argument('--kmax', type=int, default=1, help='maximum number of causal variants')

args = parser.parse_args()
database_dir = args.database_dir
sumstats = args.summary_stats
submitgap = args.submit_interval
checkgap = args.check_interval
memory = args.memory
kstart = args.kstart
kmax = args.kmax

rerunlimts = 3

# load functions
def get_checkline(filepath='dir/to/file.tsv'):
	f = open(filepath)
	checkline = f.readline()
	f.close()

	return checkline

def submit_a_chromosome(chr, memory=str(2), kstart=str(1), kmax=str(5)):
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
	# it cannot overwrite .err & .out
	if os.path.exists(tempdir + jobname + '.err'): os.remove(tempdir + jobname + '.err')
	if os.path.exists(tempdir + jobname + '.out'): os.remove(tempdir + jobname + '.out')

	memory = str(memory)
	os.system("gsub -m " + memory + " -n " + jobname + " -d " + tempdir + " -q production-rh74 'python -m cProfile -o " + tempdir + chrname + ".profile postgap/POSTGAP.py --database_dir " + database_dir + " --summary_stats " + tempdir + chrfile + " --hdf5 /nfs/production/panda/ensembl/funcgen/eqtl/GTEx.V6.88_38.cis.eqtls.h5 --sqlite /nfs/production/panda/ensembl/funcgen/eqtl/GTEx.V6.88_38.cis.eqtls.h5.sqlite3 --bayesian --kstart " + kstart + " --kmax " + kmax + " --output " + tempdir + chrname + "_res.txt --output2 " + tempdir + chrname + "_op2.txt -g' -r")
	
	print(time.strftime("%H:%M:%S", time.localtime()) + ', chr' + chr + ' is submitted with ' + memory + 'G memory assigned!')

def start_a_GWAS_sum_stats(chrlist, submitgap=90, memory=str(2), kstart=str(1), kmax=str(5)):
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
		submit_a_chromosome(chr=n, memory=str(memory), kstart=str(kstart), kmax=str(kmax))
	
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
	increase_mem = False

	# if jobstatus.txt is an empty file, this job isn't running; else, running
	os.system('bjobs -J ' + jobname + ' | cat > ' + tempdir + 'jobstatus.txt')

	# if this job is not running, check whether it finished successfully
	if os.path.getsize(tempdir + 'jobstatus.txt') == 0:
		outtext = read_job_sum_file(path=tempdir, whichjob=jobname, sumtype='out')

		if 'Successfully completed' in outtext[-23]:
			success = True
		elif 'Exited with' in outtext[-23]: #'Exited with exit code' or 'Exited with signal termination'
			repeat = True

			# get Delta Memory
			if float(outtext[-15].split()[-2]) < 0:
				increase_mem = True
				print('chr ' + chr + ': memory limit exceeded, will increase memory to re-run')
		else: # generate a NaN value, which should never happen
			success = float('NaN')
			print('something should never happen happened to chr ' + chr + '! that line was "' + outtext[-23].rstrip() + '" when checking')
	
	os.remove(tempdir + 'jobstatus.txt')

	return success, repeat, increase_mem

def track_error(chr):
	"""
	tempdir, fnstart, starttime take global varables
	"""
	jobname = fnstart + '-chr' + chr + '_' + starttime
	errtext = read_job_sum_file(path=tempdir, whichjob=jobname, sumtype='err')

	if not errtext:
		print('chr ' + chr + ': empty error file')
	elif 'HTTPError' in errtext[-1]:
		print('chr ' + chr + ': http error, no need to increase memory to re-run')
	elif 'ConnectionError' in errtext[-1]:
		print('chr ' + chr + ': connection error, no need to increase memory to re-run')
	elif 'IOError' in errtext[-1]:
		print('chr' + chr + ' gave ' + errtext[-1].rstrip() + '. performance issue with the filesystem, just try again')
	elif 'KeyboardInterrupt\n' in errtext:
		print('chr' + chr + ': KeyboardInterrupt, increase_mem should be True')
	elif 'Terminated\n' in errtext:
		print('chr' + chr + ': Terminated, increase_mem should be True')
	else:
		print('warnings, warnings, warnings! chr' + chr + ' has an unknown error situation!')
		sys.exit('unknown error situation, please debug')

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

			# just wait for 1 min, in case the stdout/stderr file has not been completed
			time.sleep(60)

			sumdf.at[i, 'success'], repeat, increase_mem = check_a_chromosome_success(chr=chr)

			# when one finished successfully, record relevant parameters
			if sumdf.loc[i, 'success'] is True:
				jobname = fnstart + '-chr' + chr + '_' + starttime
				outtext = read_job_sum_file(path=tempdir, whichjob=jobname, sumtype='out')
				
				sumdf.ix[i, 'run_time':'ave_mem'] = [outtext[-11].split()[-2], outtext[-19].split()[-2], outtext[-18].split()[-2] if 'MB' in outtext[-18] else 'N/A', outtext[-17].split()[-2] if 'MB' in outtext[-17] else 'N/A']

				print('chr' + chr + ' finished! good good')

				# once succeeds, remove the chromosomes GWAS summary statistics file to save space
				chrfile = fnstart + '-chr' + chr + '.' + filetype
				os.remove(tempdir + chrfile)

			if repeat:
				sumdf.at[i, 'repeats'] += 1
				reruntimes = sumdf.loc[i, 'repeats']
				
				if reruntimes <= rerunlimts:
					track_error(chr=chr)
					print('chr' + chr + '\'s re-run round ' + str(reruntimes) + ' starts')
					submit_a_chromosome(chr=chr, memory=str(int(init_mem) * (1 + (reruntimes if increase_mem else 0))), kstart=str(kstart), kmax=str(kmax))

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

def merge_results(chrlist):
	"""
	filedir, tempdir, fnstart, starttime take global varables
	"""
	resfile = filedir + fnstart + '_results_' + starttime + '.txt'

	# results file could have only column names
	for n in chrlist:
		chrres = tempdir + fnstart + '-chr' + n + '_res.txt'

		if n == chrlist[0]:
			os.system('cat ' + chrres + ' > ' + resfile)
		else:
			os.system('tail -n+2 ' + chrres + ' >> ' + resfile)

def merge_output2(chrlist):
	"""
	filedir, tempdir, fnstart, starttime take global varables
	"""
	op2file = filedir + fnstart + '_output2_' + starttime + '.txt'

	# output2 file could be an empty file, and has no header
	os.system('echo -e Gene_ID$"\t"Cluster_description$"\t"SNP_ID$"\t"CLPP_at_that_SNP$"\t"Tissue$"\t"CLPP_over_the_whole_cluster > ' + op2file)
	
	for n in chrlist:
		chrop2 = tempdir + fnstart + '-chr' + n + '_op2.txt'

		os.system('cat ' + chrop2 + ' >> ' + op2file)
		os.system('tail -c1 < ' + op2file + ' | read -r _ || echo ''>> ' + op2file)

#0. preparations
# get basic info about the GWAS summary statistics
filename = sumstats.split('/')[-1]
filedir = sumstats.replace(filename, '') if len(sumstats.split('/')) > 1 else './'
filetype = filename.split('.')[-1] if len(filename.split('.')) > 1 else '' #eg:tsv
fnstart = filename.replace('.' + filetype, '') #eg:WojcikGL
print(filename + ' is in the ' + ('current directory' if filedir == './' else ('directory ' + filedir)) + '\nstart pre-checking for ' + fnstart)

# create a temporary folder for all the intermediate files that generated during the analysis
tempdir = filedir + fnstart + '_' + starttime + '/'
os.makedirs(tempdir)

#1. check the filetype of the GWAS summary statistics
#a. if it is an excel file, convert it into a tsv file
if filetype in ['xls', 'xlsx', 'xlsm', 'xlsb']:
	pd.read_excel(sumstats).to_csv(tempdir + fnstart + '.tsv', header=True, index=False, sep='\t')
	print(filename + ' is an excel file, convert it into a tab-separated file, which is saved in dir ' + tempdir)
	# update filetype and filename
	filetype = 'tsv'
	filename = fnstart + '.' + filetype
elif filetype == 'out':
	#b. if it is ".out", change it to ".tsv". Note: conflict with the ".out" file from job summary
	os.rename(sumstats, filedir + fnstart + '.tsv')
	print('change filename from ' + filename + ' in ' + filedir + ' to ' + fnstart + '.tsv')
	# update filetype and filename
	filetype = 'tsv'
	filename = fnstart + '.' + filetype

	shutil.copyfile(filedir + filename, tempdir + filename)
else:
	if filetype == '':
		filetype = 'tsv'
		filename = fnstart + '.' + filetype

	shutil.copyfile(sumstats, tempdir + filename)

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
			os.remove(tempdir + fnstart + '_backup.' + filetype)

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

start_a_GWAS_sum_stats(chrlist=chrlist, submitgap=submitgap, memory=str(memory), kstart=str(kstart), kmax=str(kmax))

#5. keep tracing until all the jobs completed successfully, and make a record of running/CPU time and max/average memory
time.sleep(300)
print('start checking whether all the jobs completed successfully ...')

sumdf = trace_chromosomes_status(chrlist=chrlist, init_mem=memory, checkgap=checkgap, rerunlimts=rerunlimts)
sumdf.to_csv(filedir + fnstart + '_jobsum_' + starttime + '.txt', header=True, index=False, sep='\t', na_rep='N/A')

#6. combine all the results into one file
print('wrap all the results together, and clean up')

merge_results(chrlist=chrlist)
merge_output2(chrlist=chrlist)

# remove all the intermediate files
shutil.rmtree(tempdir)
