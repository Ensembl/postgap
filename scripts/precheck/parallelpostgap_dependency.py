# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a script file containing all the functions to run postgap.py in parallel:
* default: split by chromosomes
* optional: split by clusters
"""

import os
import pandas as pd
import shutil
import sys
import time

def get_basic_info(sumstats, starttime):
	"""
	get basic info of a GWAS summary statistics, and copy/convert it to a temporary directory
	"""
	gwasfile = sumstats.split('/')[-1]
	filedir = sumstats.replace(gwasfile, '') if len(sumstats.split('/')) > 1 else './'
	filetype = gwasfile.split('.')[-1] if len(gwasfile.split('.')) > 1 else ''
	fnstart = gwasfile.replace('.' + filetype, '')
	print('original file ' + gwasfile + ' is in the ' + ('current directory' if filedir == './' else ('directory ' + filedir)) + '\nstart pre-checking for ' + fnstart)

	# create a temporary folder for all the intermediate files that generated during the analysis
	tempdir = filedir + fnstart + '_' + starttime + '/'
	if not os.path.exists(tempdir):
		os.makedirs(tempdir)
		print('temporary directory ' + tempdir + ' is created to save intermediate files')

	# check the filetype of the GWAS summary statistics
		#a. if it is an excel file, convert it into a tsv file
	if filetype in ['xls', 'xlsx', 'xlsm', 'xlsb']:
		pd.read_excel(sumstats).to_csv(tempdir + fnstart + '.tsv', header=True, index=False, sep='\t')
		print(gwasfile + ' is an excel file, convert it to a tab-separated file and save the new file to ' + tempdir)
		# update filetype and gwasfile
		filetype = 'tsv'
		gwasfile = fnstart + '.' + filetype
	elif filetype == 'out':
		#b. if the suffix is ".out", change it to ".tsv". Note: conflict with the ".out" file from job summary
		os.rename(sumstats, filedir + fnstart + '.tsv')
		print('change filename from ' + gwasfile + ' in ' + filedir + ' to ' + fnstart + '.tsv')
		# update filetype and gwasfile
		filetype = 'tsv'
		gwasfile = fnstart + '.' + filetype

		shutil.copyfile(filedir + gwasfile, tempdir + gwasfile)
	else:
		#c. if it has no suffix, add '.tsv'
		if filetype == '':
			filetype = 'tsv'
			gwasfile = fnstart + '.' + filetype

		shutil.copyfile(sumstats, tempdir + gwasfile)

	print('new file ' + gwasfile + ' is in temporary directory ' + tempdir)
	return gwasfile, fnstart, filetype, tempdir

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

def get_clusterlist(chrlist, tempdir, fnstart):
	"""
	create a list of all the clusters
	"""
	clusterlist = []
	for chr in chrlist:
		clusterlist += os.listdir(tempdir + fnstart + '-chr' + chr + '_clusters/')
	
	return clusterlist

def submit_a_chromosome(database_dir, tempdir, chrfile, memory, kstart, kmax, split_cluster):
	"""
	submit a job to run POSTGAP for a chromosome, either to split by clusters or not
	Note:
	* type(memory/kstart/kmax): 'str'
	* split_cluster: True / False
	* chrfile = fnstart + '-chr' + chr + '.' + filetype
	* chrname = fnstart + '-chr' + chr
	"""
	starttime = tempdir.split('_')[-1].replace('/', '')
	fnstart = tempdir.split('/')[-2].replace('_' + starttime, '')
	filetype = chrfile.split('.')[-1]
	chrname = chrfile.replace('.' + filetype, '')
	job_chr = chrname + ('_spcl_' if split_cluster else '_') + starttime

	# res & op2 may exist when re-run jobs
	if os.path.exists(tempdir + chrname + '_res.txt'): os.remove(tempdir + chrname + '_res.txt')
	if os.path.exists(tempdir + chrname + '_op2.txt'): os.remove(tempdir + chrname + '_op2.txt')
	# it cannot overwrite .err & .out
	if os.path.exists(tempdir + job_chr + '.err'): os.remove(tempdir + job_chr + '.err')
	if os.path.exists(tempdir + job_chr + '.out'): os.remove(tempdir + job_chr + '.out')

	eqtl_db_part = " --hdf5 /nfs/production/panda/ensembl/funcgen/eqtl/GTEx.V6.88_38.cis.eqtls.h5 --sqlite /nfs/production/panda/ensembl/funcgen/eqtl/GTEx.V6.88_38.cis.eqtls.h5.sqlite3"
		
	if split_cluster:
		# if split by clusters, only a folder needed to save intermediate cluster files, no results or output2 will generate now
		cluster_dir = tempdir + chrname + '_clusters/'

		if not os.path.exists(cluster_dir): os.makedirs(cluster_dir)
		
		outfile_part = " --cluster_dir " + cluster_dir
	else:
		outfile_part = " --output " + tempdir + chrname + "_res.txt --output2 " + tempdir + chrname + "_op2.txt"

	POSTGAP_setting = "'python postgap/POSTGAP.py --database_dir " + database_dir + " --summary_stats " + tempdir + chrfile + eqtl_db_part + " --bayesian --kstart " + kstart + " --kmax " + kmax + outfile_part + " -g'"
	
	memory = str(memory)
	os.system("gsub -m " + memory + " -n " + job_chr + " -d " + tempdir + " -q production-rh74 " + POSTGAP_setting + " -r")
	
	print(time.strftime("%H:%M:%S", time.localtime()) + ', ' + chrname + ' is submitted with ' + memory + 'G memory assigned!')

def start_a_GWAS_sum_stats(database_dir, tempdir, gwasfile, chrlist, submitgap, memory, kstart, kmax, split_cluster):
	"""
	start to analyse a GWAS summary statistics:
	* split the original GWAS summary statistics by chromosomes
	* submit a job to run the pipeline for each, no matter splitting by clusters or not
	Note:
	* type(submitgap): 'int'
	* type(memory/kstart/kmax): 'str'
	* split_cluster: True / False
	"""
	filetype = gwasfile.split('.')[-1] # now there's a suffix for sure
	fnstart = gwasfile.replace('.' + filetype, '')

	# get which colnumn is variant_id
	checkline = get_firstline(tempdir + gwasfile)
	ncol_rsid = str(checkline.split().index('variant_id') + 1)
	print('column ' + ncol_rsid + ' is SNP rsid, start splitting by chromosomes and running pipeline' + (' with splitting by clusters' if split_cluster else '') + ':')

	for n in chrlist:
		# add a sleeping time between jobs
		if n is not chrlist[0]: time.sleep(submitgap)

		# create a new file for each chromosome, and run the analysis
		chrfile = fnstart + '-chr' + n + '.' + filetype
		os.system("(head -n1 " + tempdir + gwasfile + " && awk 'NR==FNR {a[$3]; next} $" + ncol_rsid + " in a {print $0}' yalan/variants/variants_rsID_coords-chr" + n + ".tsv " + tempdir + gwasfile + ") > " + tempdir + chrfile)
		submit_a_chromosome(database_dir, tempdir, chrfile, str(memory), str(kstart), str(kmax), split_cluster)
	
	# remove the original GWAS summary statistics file in the tempdir to save space
	os.remove(tempdir + gwasfile)

def read_job_sum_file(path='dir/', jobname='jobname', sumtype='out'):
	"""
	load the stdout/stderr file of a job
	"""
	f = open(path + jobname + '.' + sumtype)
	text = f.readlines()
	f.close()

	return text

def check_job_success(path, jobname):
	"""
	check whether a job is finished successfully
	"""
	success = False
	repeat = False
	increase_mem = False

	# if jobstatus.txt is an empty file, this job isn't running; otherwise, running
	os.system('bjobs -J ' + jobname + ' | cat > ' + path + 'jobstatus.txt')

	# if this job is not running, check whether it finished successfully
	if os.path.getsize(path + 'jobstatus.txt') == 0:
		outtext = read_job_sum_file(path, jobname, 'out')

		if 'Successfully completed' in outtext[-23]:
			success = True
		elif 'Exited with' in outtext[-23]: #'Exited with exit code' or 'Exited with signal termination'
			repeat = True

			# decide whether to increase memory or not: Delta Memory < 0 / Total Requested Memory < Max Memory
			try:
				increase_mem = (float(outtext[-15].split()[-2]) < 0 or float(outtext[-16].split()[-2]) < float(outtext[-18].split()[-2]))
				if increase_mem: print(jobname + ': memory limit exceeded, will increase memory to re-run')
			except:
				print('Oops!', sys.exc_info()[0], 'occurred when reading "' + jobname + '.out", and when checking:\n' + outtext[-15] + outtext[-16] + outtext[-16] + '\twill make a copy, this is rare')
				shutil.copyfile(path + jobname + '.out', path + jobname + '_backup.out')

				increase_mem = True
				print(jobname + ': something weird, just increase memory to re-run in any case')
		else: # something weird occurs, which should never happen
			repeat = True
			print('something should never happen happened to ' + jobname + '! that line was "' + outtext[-23].rstrip() + '" when checking')
	
	os.remove(path + 'jobstatus.txt')

	return success, repeat, increase_mem

def track_job_error(path, jobname):
	"""
	specify the type of error that failed the analysis
	"""
	errtext = read_job_sum_file(path, jobname, 'err')

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

def submit_a_cluster(database_dir, cluster_dir, cluster, chrfile, memory, kstart, kmax):
	"""
	submit a job to run POSTGAP for a cluster
	Note:
	* type(memory/kstart/kmax): 'str'
	* cluster_dir = tempdir + chrname + '_clusters/'
	"""
	tempdir = '/'.join(cluster_dir.split('/')[:-2]) + '/'
	starttime = tempdir.split('_')[-1].replace('/', '')
	fnstart = tempdir.split('/')[-2].replace('_' + starttime, '')

	clustername = fnstart + '-' + cluster
	clusterfile = cluster_dir + cluster
	job_cluster = fnstart + '-' + cluster + '_' + starttime
	
	# res & op2 may exist when re-run jobs
	if os.path.exists(tempdir + clustername + '_res.txt'): os.remove(tempdir + clustername + '_res.txt')
	if os.path.exists(tempdir + clustername + '_op2.txt'): os.remove(tempdir + clustername + '_op2.txt')
	# it cannot overwrite .err & .out
	if os.path.exists(tempdir + job_cluster + '.err'): os.remove(tempdir + job_cluster + '.err')
	if os.path.exists(tempdir + job_cluster + '.out'): os.remove(tempdir + job_cluster + '.out')

	eqtl_db_part = " --hdf5 /nfs/production/panda/ensembl/funcgen/eqtl/GTEx.V6.88_38.cis.eqtls.h5 --sqlite /nfs/production/panda/ensembl/funcgen/eqtl/GTEx.V6.88_38.cis.eqtls.h5.sqlite3"
	outfile_part = " --output " + tempdir + clustername + "_res.txt --output2 " + tempdir + clustername + "_op2.txt"

	POSTGAP_setting = "'python postgap/POSTGAP.py --database_dir " + database_dir + " --cluster_file " + clusterfile + " --summary_stats " + tempdir + chrfile + eqtl_db_part + " --bayesian --kstart " + kstart + " --kmax " + kmax + outfile_part + " -g'"
	
	memory = str(memory)
	os.system("gsub -m " + memory + " -n " + job_cluster + " -d " + tempdir + " -q production-rh74 " + POSTGAP_setting + " -r")
	
	print(time.strftime("%H:%M:%S", time.localtime()) + ', ' + cluster + ' is submitted with ' + memory + 'G memory assigned!')

def continue_chromosome_clusters(database_dir, tempdir, chrfile, submitgap, memory, kstart, kmax):
	"""
	start to analyse all the clusters of a chromosome:
	* submit a job to run the pipeline for each
	Note:
	* type(submitgap): 'int'
	* type(memory/kstart/kmax): 'str'
	"""
	chr = chrfile.split('-chr')[-1].split('.')[0]
	filetype = chrfile.split('.')[-1]
	fnstart = chrfile.split('-')[0]

	cluster_dir = tempdir + fnstart + '-chr' + chr + '_clusters/'
	chrclusters = os.listdir(cluster_dir)

	for cluster in chrclusters:
		if cluster is not chrclusters[0]: time.sleep(submitgap)

		submit_a_cluster(database_dir, cluster_dir, cluster, chrfile, str(memory), str(kstart), str(kmax))

	print('chr' + chr + ' finished submitting all the cluster jobs')

def trace_chromosomes_status(database_dir, tempdir, filetype, chrlist, submitgap, checkgap, rerunlimts, init_mem, kstart, kmax, split_cluster):
	"""
	keep checking and rerunning (if necessary) until all the chromosome jobs are finished successfully, and submit cluster jobs if split by clusters
	Note:
	* type(submitgap, checkgap, rerunlimts, init_mem): 'int'
	* type(kstart/kmax): 'str'
	* split_cluster: True / False
	"""
	starttime = tempdir.split('_')[-1].replace('/', '')
	fnstart = tempdir.split('/')[-2].replace('_' + starttime, '')

	sumdf = pd.DataFrame(columns=['chr', 'success', 'repeats', 'run_time', 'cpu_time', 'max_mem', 'ave_mem'], index=range(len(chrlist)))
	
	sumdf['chr'] = chrlist
	sumdf['repeats'] = [0] * len(chrlist)
	
	# keep checking all the unfinished jobs
	while not all(sumdf['success'] == True):
		for i,row in sumdf.loc[sumdf['success'] != True].iterrows():
			chr = row['chr']
			chrfile = fnstart + '-chr' + chr + '.' + filetype
			job_chr = fnstart + '-chr' + chr + ('_spcl_' if split_cluster else '_') + starttime

			# just wait for 1 min, in case the stdout/stderr file has not been completed
			time.sleep(60)
			sumdf.at[i, 'success'], repeat, increase_mem = check_job_success(tempdir, job_chr)

			# when one finished successfully, record relevant parameters
			if sumdf.loc[i, 'success'] is True:
				outtext = read_job_sum_file(path=tempdir, jobname=job_chr, sumtype='out')
				
				sumdf.ix[i, 'run_time':'ave_mem'] = [outtext[-11].split()[-2], outtext[-19].split()[-2], outtext[-18].split()[-2] if 'MB' in outtext[-18] else 'N/A', outtext[-17].split()[-2] if 'MB' in outtext[-17] else 'N/A']

				if split_cluster:
					#once splitting by clusters is done, submit a job to run the pipeline for each cluster
					print('chr' + chr + ' completed creating cluster files, and start to run POSTGAP for clusters')
					continue_chromosome_clusters(database_dir, tempdir, chrfile, submitgap, str(init_mem), str(kstart), str(kmax))
				else:
					print('chr' + chr + ' finished! good good')

			if repeat:
				sumdf.at[i, 'repeats'] += 1
				reruntimes = sumdf.loc[i, 'repeats']
				
				if reruntimes <= rerunlimts:
					print('chr' + chr + ' gave: ' + track_job_error(tempdir, job_chr))

					print('chr' + chr + '\'s re-run round ' + str(reruntimes) + ' starts')
					submit_a_chromosome(database_dir, tempdir, chrfile, str(init_mem * (1 + (reruntimes if increase_mem else 0))), str(kstart), str(kmax), split_cluster)

		if all(sumdf['success'] == True):
			print('congratulations!!! all the chr' + ('_spcl' if split_cluster else '') + ' jobs finished successfully!')
			break
		elif any((sumdf['success'] == False) & (sumdf['repeats'] > rerunlimts)):
			print('/'.join(('chr' + n) for n in sumdf.loc[(sumdf['success'] == False) & (sumdf['repeats'] > rerunlimts), 'chr']) + ' has already re-ran for ' + str(rerunlimts) + ' times, some consistent error occurs!')
			sys.exit('consistent error occurs, please debug')

		time.sleep(checkgap)
	
	return sumdf

def trace_clusters_status(database_dir, tempdir, filetype, clusterlist, checkgap, rerunlimts, init_mem, kstart, kmax):
	"""
	keep checking and rerunning (if necessary) until all the cluster jobs are finished successfully
	Note:
	* type(checkgap, rerunlimts, init_mem): 'int'
	* type(kstart/kmax): 'str'
	"""
	starttime = tempdir.split('_')[-1].replace('/', '')
	fnstart = tempdir.split('/')[-2].replace('_' + starttime, '')

	sumdf = pd.DataFrame(columns=['cluster', 'success', 'repeats', 'fail_fm', 'run_time', 'cpu_time', 'max_mem', 'ave_mem'], index=range(len(clusterlist)))
	sumdf['cluster'] = clusterlist
	sumdf['repeats'] = [0] * len(clusterlist)
	
	# keep checking all the unfinished jobs
	while not all(sumdf['success'] == True):
		for i,row in sumdf.loc[sumdf['success'] != True].iterrows():
			cluster = row['cluster'] #'GWAS_Cluster_21:44931821-45100521'
			job_cluster = fnstart + '-' + cluster + '_' + starttime

			chr = cluster.split(':')[0].split('_')[-1]
			cluster_dir = tempdir + fnstart + '-chr' + chr + '_clusters/'
			chrfile = fnstart + '-chr' + chr + '.' + filetype

			# just wait for 1 min, in case the stdout/stderr file has not been completed
			time.sleep(60)
			sumdf.at[i, 'success'], repeat, increase_mem = check_job_success(tempdir, job_cluster)

			# when one job finished successfully, record relevant parameters
			if sumdf.loc[i, 'success'] is True:
				outtext = read_job_sum_file(path=tempdir, jobname=job_cluster, sumtype='out')
				
				fail_fm = any(['failed finemapping\n' in x for x in outtext])
				sumdf.ix[i, 'fail_fm':'ave_mem'] = [fail_fm, outtext[-11].split()[-2], outtext[-19].split()[-2], outtext[-18].split()[-2] if 'MB' in outtext[-18] else 'N/A', outtext[-17].split()[-2] if 'MB' in outtext[-17] else 'N/A']

				print(cluster + ' finished! ' + ('but failed finemapping' if fail_fm else 'good good'))

			if repeat:
				last_mem = int(float(outtext[-16].split()[-2]) / 1000) # 5000.00 MB -> 5 GB

				sumdf.at[i, 'repeats'] = last_mem / init_mem
				reruntimes = sumdf.loc[i, 'repeats']
				
				if reruntimes <= rerunlimts:
					print(cluster + ' gave: ' + track_job_error(tempdir, job_cluster))

					print(cluster + '\'s re-run round ' + str(reruntimes) + ' starts')
					submit_a_cluster(database_dir, cluster_dir, cluster, chrfile, str(last_mem + (init_mem if increase_mem else 0)), str(kstart), str(kmax))

		if all(sumdf['success'] == True):
			print('congratulations!!! all the cluster jobs finished successfully!')
			break
		elif any((sumdf['success'] == False) & (sumdf['repeats'] > rerunlimts)):
			print('/'.join(c for c in sumdf.loc[(sumdf['success'] == False) & (sumdf['repeats'] > rerunlimts), 'cluster']) + ' has already re-ran for ' + str(rerunlimts) + ' times, some consistent error occurs!')
			sys.exit('consistent error occurs, please debug')

		time.sleep(checkgap)
	
	return sumdf

def merge_results(objectlist, tempdir, kstart, kmax, split_cluster):
	"""
	merge all the results files into a big one
	Note:
	* type(kstart/kmax): 'str'
	* split_cluster: True / False
	"""
	starttime = tempdir.split('_')[-1].replace('/', '')
	fnstart = tempdir.split('/')[-2].replace('_' + starttime, '')

	resfile = '/'.join(tempdir.split('/')[:-2]) + '/' + fnstart + '_results_ks' + str(kstart) + '_km' + str(kmax) + ('_spcl_' if split_cluster else '_') + starttime + '.txt'

	# results file could have only column names
	for n in objectlist:
		eachres = tempdir + fnstart + '-' + (n if split_cluster else ('chr' + n)) + '_res.txt'

		if n == objectlist[0]:
			os.system('cat ' + eachres + ' > ' + resfile)
		else:
			os.system('tail -n+2 ' + eachres + ' >> ' + resfile)

def merge_output2(objectlist, tempdir, kstart, kmax, split_cluster):
	"""
	merge all the output2 files into a big one
	Note:
	* type(kstart/kmax): 'str'
	* split_cluster: True / False
	"""
	starttime = tempdir.split('_')[-1].replace('/', '')
	fnstart = tempdir.split('/')[-2].replace('_' + starttime, '')

	op2file = '/'.join(tempdir.split('/')[:-2]) + '/' + fnstart + '_output2_ks' + str(kstart) + '_km' + str(kmax) + ('_spcl_' if split_cluster else '_') + starttime + '.txt'

	# output2 file could be an empty file, and has no header
	os.system('echo -e Gene_ID$"\t"Cluster_description$"\t"SNP_ID$"\t"CLPP_at_that_SNP$"\t"Tissue$"\t"CLPP_over_the_whole_cluster > ' + op2file)
	
	for n in objectlist:
		eachop2 = tempdir + fnstart + '-' + (n if split_cluster else ('chr' + n)) + '_op2.txt'

		os.system('cat ' + eachop2 + ' >> ' + op2file)
		# add a new line character, as op2 file does not include
		os.system('tail -c1 < ' + op2file + ' | read -r _ || echo ''>> ' + op2file)
