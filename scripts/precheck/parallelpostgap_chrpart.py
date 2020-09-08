# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a script file to run postgap.py in parallel:
* default: split by chromosomes
* optional: split by clusters

Note: temporarily run the script under the home directory!!!
"""

import argparse
import os
import pandas as pd
import shutil
import sys
import time
from parallelpostgap_dependency import *

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

starttime = time.strftime("%y%m%d%H%M%S", time.localtime())
rerunlimts = 4

#1. get basic info about the GWAS summary statistics, and copy the GWAS file into the temporary directory
gwasfile, fnstart, filetype, tempdir = get_basic_info(sumstats, starttime)

#2. check if the delimiter is tab, and convert if it is one of the common ones
check_delimiter(file=tempdir + gwasfile)
# now the separator is tab for sure

#3. check three very important colnames: variant_id, beta, p-value
check_three_colnames(file=tempdir + gwasfile)

#4. split GWAS summary file by chromosomes and run postgap.py on each
chrlist = define_chrlist()
print('Note: split by clusters is ' + ('ON' if split_cluster else 'OFF'))
start_a_GWAS_sum_stats(database_dir, tempdir, gwasfile, chrlist, submitgap, str(memory), str(kstart), str(kmax), split_cluster)

#5. keep tracing until all the jobs completed successfully, and make a record of running/CPU time and max/average memory
print('start checking whether all the chromosome jobs completed successfully ...')

sumdf = trace_chromosomes_status(database_dir, tempdir, filetype, chrlist, submitgap, checkgap, rerunlimts, memory, str(kstart), str(kmax), split_cluster)
sumdf.to_csv('/'.join(tempdir.split('/')[:-2]) + '/' + fnstart + '_chrsum_ks' + str(kstart) + '_km' + str(kmax) + ('_spcl_' if split_cluster else '_') + starttime + '.txt', header=True, index=False, sep='\t', na_rep='N/A')

if split_cluster:
	print(fnstart + ' will start a new job to check all the cluster jobs')
	# when reach here, it means all the chr jobs are done, and all the cluster jobs are submitted
	os.system("gsub -m 10 -n " + fnstart + "_" + str(kstart) + str(kmax) + "_clpt -d yalan/NNRtest/ -q production-rh74 'python parallelpostgap_clusterpart.py --database_dir " + database_dir + " --check_interval " + str(checkgap) + " --memory " + str(memory) + " --kstart " + str(kstart) + " --kmax " + str(kmax) + " --tempdir " + tempdir + " --filetype " + filetype + " --rerunlimts " + str(rerunlimts) + "' -r")
else:
	#6. combine all the results into one file
	print('wrap all the results together, and clean up')

	merge_results(chrlist, tempdir, str(kstart), str(kmax), split_cluster)
	merge_output2(chrlist, tempdir, str(kstart), str(kmax), split_cluster)

	print(fnstart + ' is done!')
