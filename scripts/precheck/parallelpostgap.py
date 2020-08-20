# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a script file for preworks before running postgap.py in parallel:
* default: split by chromosomes
* optional: split by clusters
"""

import argparse
import os
import pandas as pd
import shutil
import sys
import time

# Note: temporarily run this code under the home directory!!!
starttime = time.strftime("%y%m%d%H%M%S", time.localtime())
rerunlimts = 4

# get command line arguments
parser = argparse.ArgumentParser(description='to read command line arguments')
parser.add_argument('--database_dir', type=str, default='databases/', help='directory where data files are stored')
parser.add_argument('--summary_stats', required=True, type=str, help='location of input summary statistics file')
parser.add_argument('--submit_interval', type=int, default=90, help='sleeping time (s) added to start the analysis of the next chromosome')
parser.add_argument('--check_interval', type=int, default=300, help='sleeping time (s) added to check again if all the analyses of chromosomes have been done')
parser.add_argument('--memory', type=int, default=2, help='memory of a node (GB) requested to run a job')
parser.add_argument('--kstart', type=int, default=1, help='how many causal variants to start with in the full exploration of sets')
parser.add_argument('--kmax', type=int, default=1, help='maximum number of causal variants')
parser.add_argument('--split_cluster', action = 'store_true', help='split by clusters to speed up the computation')

args = parser.parse_args()
database_dir = args.database_dir
sumstats = args.summary_stats
submitgap = args.submit_interval
checkgap = args.check_interval
memory = args.memory
kstart = args.kstart
kmax = args.kmax
split_cluster = args.split_cluster

# load functions
def get_basic_info(sumstats):
	"""
	get basic info of a GWAS summary statistics, and copy/convert it to a temporary directory
	"""
	filename = sumstats.split('/')[-1]
	filedir = sumstats.replace(filename, '') if len(sumstats.split('/')) > 1 else './'
	filetype = filename.split('.')[-1] if len(filename.split('.')) > 1 else ''
	fnstart = filename.replace('.' + filetype, '')
	print('original file ' + filename + ' is in the ' + ('current directory' if filedir == './' else ('directory ' + filedir)) + '\nstart pre-checking for ' + fnstart)

	# create a temporary folder for all the intermediate files that generated during the analysis
	tempdir = filedir + fnstart + '_' + starttime + '/'
	if not os.path.exists(tempdir):
		os.makedirs(tempdir)
		print('temporary directory ' + tempdir + ' is created to save intermediate files')

	# check the filetype of the GWAS summary statistics
		#a. if it is an excel file, convert it into a tsv file
	if filetype in ['xls', 'xlsx', 'xlsm', 'xlsb']:
		pd.read_excel(sumstats).to_csv(tempdir + fnstart + '.tsv', header=True, index=False, sep='\t')
		print(filename + ' is an excel file, convert it to a tab-separated file and save the new file to ' + tempdir)
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
		#c. if it has no suffix, add '.tsv'
		if filetype == '':
			filetype = 'tsv'
			filename = fnstart + '.' + filetype

		shutil.copyfile(sumstats, tempdir + filename)

	print('new file ' + filename + ' is in temporary directory ' + tempdir)
	return filedir, filename, fnstart, filetype, tempdir

def get_firstline(file='dir/to/file.tsv'):
	"""
	get the first line in a file
	"""
	f = open(file)
	firstline = f.readline()
	f.close()

	return firstline

def check_delimiter(file='dir/to/file.tsv'):
	"""
	check whether the delimiter of a file is tab or not, and convert it to a tab-sep file if its delimiter is one of common ones
	"""
	checkline = get_firstline(file)

	if len(checkline.split('\t')) > 1:
		print('delimiter checking passed')
	else:
		seplist=[' ', ',', ';', ':', '|']

		t = 0
		for s in seplist:
			if len(checkline.split(s)) > 1:
				print('delimiter is ' + ('space' if s == ' ' else s) + ', and change it to tab!')

				file_type = file.split('.')[-1] if len(file.split('.')) > 1 else ''
				file_rest = file.replace('.' + file_type, '')
				# replace delimiter with tab, and save into a new file
				os.rename(file, file_rest + '_backup.' + file_type)
				os.system('cat ' + file_rest + '_backup.' + file_type + ' | tr "' + s + '" "\t" > ' + file)
				os.remove(file_rest + '_backup.' + file_type)

				break
			else:
				t += 1
		if t == len(seplist):
			print('unknown delimiter, please change the delimiter to tab')
			sys.exit('error with delimiter')

def check_three_colnames(file='dir/to/file.tsv'):
	"""
	check whether a GWAS summary statistics uses 'variant_id/beta/p-value' as column names
	"""
	checkline = get_firstline(file)

	if all([impnames in checkline.split() for impnames in ['variant_id', 'beta', 'p-value']]):
		print('colnames checking passed')
	else:
		print('please make sure all the following information is included in your data, and please specify them with assigned column names:\nrsid -> "variant_id"\nbeta -> "beta"\np -> "p-value"')
		sys.exit('error with column names')

def define_chrlist():
	"""
	create a list of chromosomes
	"""
	chrlist = [str(n) for n in range(1, 23)]
	chrlist.append('O')

	return chrlist

def submit_a_chromosome(chr, memory=str(2), kstart=str(1), kmax=str(5), split_cluster=False):
	"""
	submit a job to run POSTGAP for a chromosome, either to split by clusters or not
	Note: type(memory/kstart/kmax) is 'str'
	database_dir, tempdir, fnstart, filetype, starttime take global variables
	"""
	chrname = fnstart + '-chr' + chr
	chrfile = fnstart + '-chr' + chr + '.' + filetype
	job_chr = fnstart + '-chr' + chr + ('_spcl_' if split_cluster else '_') + starttime

	# res & op2 may exist when re-run jobs
	if os.path.exists(tempdir + chrname + '_res.txt'): os.remove(tempdir + chrname + '_res.txt')
	if os.path.exists(tempdir + chrname + '_op2.txt'): os.remove(tempdir + chrname + '_op2.txt')
	if os.path.exists(tempdir + chrname + ('_spcl' if split_cluster else '') + '.profile'): os.remove(tempdir + chrname + ('_spcl' if split_cluster else '') + '.profile')
	# it cannot overwrite .err & .out
	if os.path.exists(tempdir + job_chr + '.err'): os.remove(tempdir + job_chr + '.err')
	if os.path.exists(tempdir + job_chr + '.out'): os.remove(tempdir + job_chr + '.out')

	profile_part = " -m cProfile -o " + tempdir + chrname + ("_spcl" if split_cluster else "") + ".profile"
	eqtl_db_part = " --hdf5 /nfs/production/panda/ensembl/funcgen/eqtl/GTEx.V6.88_38.cis.eqtls.h5 --sqlite /nfs/production/panda/ensembl/funcgen/eqtl/GTEx.V6.88_38.cis.eqtls.h5.sqlite3"
		
	if split_cluster:
		# if split by clusters, only a folder needed to save intermediate cluster files, no results or output2 will generate now
		cluster_dir = tempdir + fnstart + '-chr' + chr + '_clusters/'

		if not os.path.exists(cluster_dir): os.makedirs(cluster_dir)
		
		outfile_part = " --cluster_dir " + cluster_dir
	else:
		outfile_part = " --output " + tempdir + chrname + "_res.txt --output2 " + tempdir + chrname + "_op2.txt"

	POSTGAP_setting = "'python" + profile_part + " postgap/POSTGAP.py --database_dir " + database_dir + " --summary_stats " + tempdir + chrfile + eqtl_db_part + " --bayesian --kstart " + kstart + " --kmax " + kmax + outfile_part + " -g'"
	
	memory = str(memory)
	os.system("gsub -m " + memory + " -n " + job_chr + " -d " + tempdir + " -q production-rh74 " + POSTGAP_setting + " -r")
	
	print(time.strftime("%H:%M:%S", time.localtime()) + ', chr' + chr + ' is submitted with ' + memory + 'G memory assigned!')

def start_a_GWAS_sum_stats(file, chrlist, split_cluster, submitgap=90, memory=str(2), kstart=str(1), kmax=str(5)):
	"""
	start to analyse a GWAS summary statistics:
	* split the original GWAS summary statistics by chromosomes
	* submit a job to run POSTGAP for each, either to split by clusters or not
	Note: type(memory/kstart/kmax) is 'str'
	tempdir, fnstart, filetype take global variables
	"""
	# get which colnumn is variant_id
	checkline = get_firstline(file)
	ncol_rsid = str(checkline.split().index('variant_id') + 1)
	print('column ' + ncol_rsid + ' is SNP rsid, start splitting by chromosomes and running pipeline' + (' with splitting by clusters' if split_cluster else '') + ':')

	for n in chrlist:
		# add a sleeping time between jobs
		if n is not chrlist[0]: time.sleep(submitgap)

		# create a new file for each chromosome, and run postgap.py
		chrfile = fnstart + '-chr' + n + '.' + filetype
		os.system("(head -n1 " + file + " && awk 'NR==FNR {a[$3]; next} $" + ncol_rsid + " in a {print $0}' yalan/variants/variants_rsID_coords-chr" + n + ".tsv " + file + ") > " + tempdir + chrfile)
		submit_a_chromosome(n, str(memory), str(kstart), str(kmax), split_cluster)
	
	# remove the original GWAS summary statistics file in the tempdir to save space
	os.remove(file)

def read_job_sum_file(path='dir/', jobname='jobname', sumtype='out'):
	"""
	load the stdout/stderr file of a job
	"""
	f = open(path + jobname + '.' + sumtype)
	text = f.readlines()
	f.close()

	return text

def check_job_success(jobname):
	"""
	check whether a job is finished successfully
	tempdir takes global variables
	"""
	success = False
	repeat = False
	increase_mem = False

	# if jobstatus.txt is an empty file, this job isn't running; else, running
	os.system('bjobs -J ' + jobname + ' | cat > ' + tempdir + 'jobstatus.txt')

	# if this job is not running, check whether it finished successfully
	if os.path.getsize(tempdir + 'jobstatus.txt') == 0:
		outtext = read_job_sum_file(path=tempdir, jobname=jobname, sumtype='out')

		if 'Successfully completed' in outtext[-23]:
			success = True
		elif 'Exited with' in outtext[-23]: #'Exited with exit code' or 'Exited with signal termination'
			repeat = True

			# get Delta Memory to decide whether to increase memory or not
			if float(outtext[-15].split()[-2]) < 0:
				increase_mem = True
				print(jobname + ': memory limit exceeded, will increase memory to re-run')
		else: # something weird occurs, which should never happen
			repeat = True
			#success = float('NaN') #any(sumdf['success'].isnull())
			print('something should never happen happened to ' + jobname + '! that line was "' + outtext[-23].rstrip() + '" when checking')
	
	os.remove(tempdir + 'jobstatus.txt')

	return success, repeat, increase_mem

def track_job_error(jobname):
	"""
	specify the type of error that failed the analysis
	tempdir takes global variables
	"""
	errtext = read_job_sum_file(path=tempdir, jobname=jobname, sumtype='err')

	if not errtext:
		return 'empty error file'
	elif 'HTTPError' in errtext[-1]:
		return 'http error, no need to increase memory to re-run'
	elif 'ConnectionError' in errtext[-1]:
		return 'connection error, no need to increase memory to re-run'
	elif 'IOError' in errtext[-1]:
		return errtext[-1].rstrip() + '. performance issue with the filesystem, just try again'
	elif 'KeyboardInterrupt\n' in errtext:
		return '\'KeyboardInterrupt\', increase_mem should be True'
	elif 'Terminated\n' in errtext:
		return '\'Terminated\', increase_mem should be True'
	else:
		print('warnings, warnings, warnings! ' + jobname + ' has an unknown error situation!')
		sys.exit('unknown error situation, please debug')

def submit_a_cluster(cluster_dir, cluster, memory=str(2), kstart=str(1), kmax=str(5)):
	"""
	submit a job to run POSTGAP for a cluster
	Note: type(memory/kstart/kmax) is 'str'
	database_dir, tempdir, fnstart, filetype, starttime take global variables
	"""
	clustername = fnstart + '-' + cluster
	clusterfile = cluster_dir + cluster
	job_cluster = fnstart + '-' + cluster + '_' + starttime
	
	chrfile = fnstart + cluster_dir.split('-')[-1].replace('_clusters', '') + '.' + filetype

	# res & op2 may exist when re-run jobs
	if os.path.exists(tempdir + clustername + '_res.txt'): os.remove(tempdir + clustername + '_res.txt')
	if os.path.exists(tempdir + clustername + '_op2.txt'): os.remove(tempdir + clustername + '_op2.txt')
	if os.path.exists(tempdir + clustername + '.profile'): os.remove(tempdir + clustername + '.profile')
	# it cannot overwrite .err & .out
	if os.path.exists(tempdir + job_cluster + '.err'): os.remove(tempdir + job_cluster + '.err')
	if os.path.exists(tempdir + job_cluster + '.out'): os.remove(tempdir + job_cluster + '.out')

	profile_part = " -m cProfile -o " + tempdir + clustername + ".profile"
	eqtl_db_part = " --hdf5 /nfs/production/panda/ensembl/funcgen/eqtl/GTEx.V6.88_38.cis.eqtls.h5 --sqlite /nfs/production/panda/ensembl/funcgen/eqtl/GTEx.V6.88_38.cis.eqtls.h5.sqlite3"
	outfile_part = " --output " + tempdir + clustername + "_res.txt --output2 " + tempdir + clustername + "_op2.txt"

	POSTGAP_setting = "'python" + profile_part + " postgap/POSTGAP.py --database_dir " + database_dir + " --cluster_file " + clusterfile + " --summary_stats " + tempdir + chrfile + eqtl_db_part + " --bayesian --kstart " + kstart + " --kmax " + kmax + outfile_part + " -g'"
	
	memory = str(memory)
	os.system("gsub -m " + memory + " -n " + job_cluster + " -d " + tempdir + " -q production-rh74 " + POSTGAP_setting + " -r")
	
	print(time.strftime("%H:%M:%S", time.localtime()) + ', ' + cluster + ' is submitted with ' + memory + 'G memory assigned!')

def continue_chromosome_clusters(chr, submitgap=90, memory=str(2), kstart=str(1), kmax=str(5)):
	"""
	start to analyse all the cluster files from a chromosome:
	* submit a job to run POSTGAP for each
	Note: type(memory/kstart/kmax) is 'str'
	tempdir, fnstart take global variables
	"""
	cluster_dir = tempdir + fnstart + '-chr' + chr + '_clusters/'
	chrclusters = os.listdir(cluster_dir)

	for cluster in chrclusters:
		if cluster is not chrclusters[0]: time.sleep(submitgap)

		submit_a_cluster(cluster_dir, cluster, str(memory), str(kstart), str(kmax))

	print('chr' + chr + ' finished submitting all the cluster jobs')
	return chrclusters

def trace_chromosomes_status(chrlist, split_cluster, init_mem=2, checkgap=300, rerunlimts=3):
	"""
	keep checking and rerunning (if necessary) until all the chromosome jobs are finished successfully, and submit cluster jobs if split by clusters
	Note: type(init_mem) is 'int'
	tempdir, fnstart, starttime take global variables
	"""
	sumdf = pd.DataFrame(columns=['chr', 'success', 'repeats', 'run_time', 'cpu_time', 'max_mem', 'ave_mem'], index=range(len(chrlist)))
	
	sumdf['chr'] = chrlist
	sumdf['repeats'] = [0] * len(chrlist)
	
	if split_cluster: allclusters = []

	while not all(sumdf['success'] == True):
		# check all unfinished jobs
		for i,row in sumdf.loc[sumdf['success'] != True].iterrows():
			chr = row['chr']
			job_chr = fnstart + '-chr' + chr + ('_spcl_' if split_cluster else '_') + starttime

			# just wait for 1 min, in case the stdout/stderr file has not been completed
			time.sleep(60)
			sumdf.at[i, 'success'], repeat, increase_mem = check_job_success(job_chr)

			# when one finished successfully, record relevant parameters
			if sumdf.loc[i, 'success'] is True:
				outtext = read_job_sum_file(path=tempdir, jobname=job_chr, sumtype='out')
				
				sumdf.ix[i, 'run_time':'ave_mem'] = [outtext[-11].split()[-2], outtext[-19].split()[-2], outtext[-18].split()[-2] if 'MB' in outtext[-18] else 'N/A', outtext[-17].split()[-2] if 'MB' in outtext[-17] else 'N/A']

				if split_cluster:
					#track if the first step splitting cluster finishes or not. once done, postgap for each cluster
					print('chr' + chr + ' completed creating cluster files, and start to run POSTGAP for clusters')
					
					chrclusters = continue_chromosome_clusters(chr, submitgap, str(memory), str(kstart), str(kmax))
					allclusters += chrclusters
				else:
					print('chr' + chr + ' finished! good good')

			if repeat:
				sumdf.at[i, 'repeats'] += 1
				reruntimes = sumdf.loc[i, 'repeats']
				
				if reruntimes <= rerunlimts:
					print('chr' + chr + ' gave: ' + track_job_error(job_chr))

					print('chr' + chr + '\'s re-run round ' + str(reruntimes) + ' starts')
					submit_a_chromosome(chr, str(int(init_mem) * (1 + (reruntimes if increase_mem else 0))), str(kstart), str(kmax), split_cluster)

		if all(sumdf['success'] == True):
			print('congratulations!!! all the chr' + ('_spcl' if split_cluster else '') + ' jobs finished successfully!')
			break
		elif any((sumdf['success'] == False) & (sumdf['repeats'] > rerunlimts)):
			print('/'.join(('chr' + n) for n in sumdf.loc[(sumdf['success'] == False) & (sumdf['repeats'] > rerunlimts), 'chr']) + ' has already re-ran for ' + str(rerunlimts) + ' times, some consistent error occurs!')
			sys.exit('consistent error occurs, please debug')

		if int(starttime[-4:]) - checkgap / 3 * 5 < int(time.strftime("%M%S", time.localtime())) < int(starttime[-4:]) + checkgap / 3 * 5:
			print('/'.join(('chr' + n) for n in sumdf.loc[sumdf['success'] != True, 'chr']) + ' unfinished, check every ' + str(checkgap) + 's')
		time.sleep(checkgap)
	
	return (sumdf, allclusters) if split_cluster else sumdf

def trace_clusters_status(clusterlist, init_mem=2, checkgap=300, rerunlimts=3):
	"""
	keep checking and rerunning (if necessary) until all the cluster jobs are finished successfully
	Note: type(init_mem) is 'int'
	tempdir, fnstart, starttime take global variables
	"""
	sumdf = pd.DataFrame(columns=['cluster', 'success', 'repeats', 'run_time', 'cpu_time', 'max_mem', 'ave_mem'], index=range(len(clusterlist)))
	sumdf['cluster'] = clusterlist
	sumdf['repeats'] = [0] * len(clusterlist)
	
	while not all(sumdf['success'] == True):
		# check all unfinished jobs
		for i,row in sumdf.loc[sumdf['success'] != True].iterrows():
			cluster = row['cluster'] #'GWAS_Cluster_21:44931821-45100521'
			job_cluster = fnstart + '-' + cluster + '_' + starttime

			chr = cluster.split(':')[0].split('_')[-1]
			cluster_dir = tempdir + fnstart + '-chr' + chr + '_clusters/'

			# just wait for 1 min, in case the stdout/stderr file has not been completed
			time.sleep(60)
			sumdf.at[i, 'success'], repeat, increase_mem = check_job_success(job_cluster)

			# when one finished successfully, record relevant parameters
			if sumdf.loc[i, 'success'] is True:
				outtext = read_job_sum_file(path=tempdir, jobname=job_cluster, sumtype='out')
				
				sumdf.ix[i, 'run_time':'ave_mem'] = [outtext[-11].split()[-2], outtext[-19].split()[-2], outtext[-18].split()[-2] if 'MB' in outtext[-18] else 'N/A', outtext[-17].split()[-2] if 'MB' in outtext[-17] else 'N/A']

				print(cluster + ' finished! good good')

			if repeat:
				sumdf.at[i, 'repeats'] += 1
				reruntimes = sumdf.loc[i, 'repeats']
				
				if reruntimes <= rerunlimts:
					print(cluster + ' gave: ' + track_job_error(job_cluster))

					print(cluster + '\'s re-run round ' + str(reruntimes) + ' starts')
					submit_a_cluster(cluster_dir, cluster, str(int(init_mem) * (1 + (reruntimes if increase_mem else 0))), str(kstart), str(kmax))

		if all(sumdf['success'] == True):
			print('congratulations!!! all the cluster jobs finished successfully!')
			break
		elif any((sumdf['success'] == False) & (sumdf['repeats'] > rerunlimts)):
			print('/'.join(c for c in sumdf.loc[(sumdf['success'] == False) & (sumdf['repeats'] > rerunlimts), 'cluster']) + ' has already re-ran for ' + str(rerunlimts) + ' times, some consistent error occurs!')
			sys.exit('consistent error occurs, please debug')

		if int(starttime[-4:]) - checkgap / 3 * 5 < int(time.strftime("%M%S", time.localtime())) < int(starttime[-4:]) + checkgap / 3 * 5:
			print('/'.join(c for c in sumdf.loc[sumdf['success'] != True, 'cluster']) + ' unfinished, check every ' + str(checkgap) + 's')
		time.sleep(checkgap)
	
	return sumdf

def merge_results(objectlist, split_cluster):
	"""
	merge all the results files into a big one
	filedir, tempdir, fnstart, starttime take global variables
	"""
	resfile = filedir + fnstart + '_results_' + starttime + '.txt'

	# results file could have only column names
	for n in objectlist:
		eachres = tempdir + fnstart + '-' + (n if split_cluster else ('chr' + n)) + '_res.txt'

		if n == objectlist[0]:
			os.system('cat ' + eachres + ' > ' + resfile)
		else:
			os.system('tail -n+2 ' + eachres + ' >> ' + resfile)

def merge_output2(objectlist, split_cluster):
	"""
	merge all the output2 files into a big one
	filedir, tempdir, fnstart, starttime take global variables
	"""
	op2file = filedir + fnstart + '_output2_' + starttime + '.txt'

	# output2 file could be an empty file, and has no header
	os.system('echo -e Gene_ID$"\t"Cluster_description$"\t"SNP_ID$"\t"CLPP_at_that_SNP$"\t"Tissue$"\t"CLPP_over_the_whole_cluster > ' + op2file)
	
	for n in objectlist:
		eachop2 = tempdir + fnstart + '-' + (n if split_cluster else ('chr' + n)) + '_op2.txt'

		os.system('cat ' + eachop2 + ' >> ' + op2file)
		# add a new line character, as op2 file does not include
		os.system('tail -c1 < ' + op2file + ' | read -r _ || echo ''>> ' + op2file)

#1. get basic info about the GWAS summary statistics, and copy the GWAS file into the temporary directory
filedir, filename, fnstart, filetype, tempdir = get_basic_info(sumstats)

#2. check if the delimiter is tab, and convert if it is one of the common ones
check_delimiter(file=tempdir + filename)
# now the separator is tab for sure

#3. check three very important colnames: variant_id, beta, p-value
check_three_colnames(file=tempdir + filename)

#4. split GWAS summary file by chromosomes and run postgap.py on each
chrlist = define_chrlist()
print('Note: split by clusters is ' + ('ON' if split_cluster else 'OFF'))
start_a_GWAS_sum_stats(tempdir + filename, chrlist, split_cluster, submitgap, str(memory), str(kstart), str(kmax))

#5. keep tracing until all the jobs completed successfully, and make a record of running/CPU time and max/average memory
print('start checking whether all the jobs completed successfully ...')

if split_cluster:
	sum_chr, allclusters = trace_chromosomes_status(chrlist, split_cluster, memory, checkgap, rerunlimts)
	# when reach here, it means all the chr jobs are done, and all the cluster jobs are submitted
	
	sum_cluster = trace_clusters_status(allclusters, memory, checkgap, rerunlimts)

	sumdf = pd.concat([sum_chr, sum_cluster], ignore_index=True)
	# TODO: change colnames - 'chr' -> 'chr/cluster'
	sumdf.to_csv(filedir + fnstart + '_spcl_jobsum_' + starttime + '.txt', header=True, index=False, sep='\t', na_rep='N/A')
else:
	sum_chr = trace_chromosomes_status(chrlist, split_cluster, memory, checkgap, rerunlimts)
	sum_chr.to_csv(filedir + fnstart + '_jobsum_' + starttime + '.txt', header=True, index=False, sep='\t', na_rep='N/A')

#6. combine all the results into one file
print('wrap all the results together, and clean up')

objectlist = allclusters if split_cluster else chrlist
merge_results(objectlist, split_cluster)
merge_output2(objectlist, split_cluster)

# remove all the intermediate files
shutil.rmtree(tempdir)
