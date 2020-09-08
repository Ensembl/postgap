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
parser.add_argument('--check_interval', type=int, default=300, help='sleeping time (s) added to check again if all the analyses of chromosomes have been done')
parser.add_argument('--memory', type=int, default=2, help='memory of a node (GB) requested to run a job')
parser.add_argument('--kstart', type=int, default=1, help='how many causal variants to start with in the full exploration of sets')
parser.add_argument('--kmax', type=int, default=1, help='maximum number of causal variants')
parser.add_argument('--tempdir',  required=True, type=str, help='')
parser.add_argument('--filetype', required=True, type=str, help='')
parser.add_argument('--rerunlimts', type=int, default=2, help='')

args = parser.parse_args()
database_dir = args.database_dir
checkgap = args.check_interval
memory = args.memory
kstart = args.kstart
kmax = args.kmax
tempdir = args.tempdir
filetype = args.filetype
rerunlimts = args.rerunlimts

starttime = tempdir.split('_')[-1].replace('/', '')
fnstart = tempdir.split('/')[-2].replace('_' + starttime, '')
split_cluster = True

#5. keep tracing until all the jobs completed successfully, and make a record of running/CPU time and max/average memory
print('start checking whether all the cluster jobs completed successfully ...')

chrlist = define_chrlist()
clusterlist = get_clusterlist(chrlist, tempdir, fnstart)

sumdf = trace_clusters_status(database_dir, tempdir, filetype, clusterlist, checkgap, rerunlimts, memory, str(kstart), str(kmax))
sumdf.to_csv('/'.join(tempdir.split('/')[:-2]) + '/' + fnstart + '_clusum_ks' + str(kstart) + '_km' + str(kmax) + ('_spcl_' if split_cluster else '_') + starttime + '.txt', header=True, index=False, sep='\t', na_rep='N/A')

#6. combine all the results into one file
print('wrap all the results together, and clean up')

merge_results(clusterlist, tempdir, str(kstart), str(kmax), split_cluster)
merge_output2(clusterlist, tempdir, str(kstart), str(kmax), split_cluster)

print(fnstart + ' is done!')
