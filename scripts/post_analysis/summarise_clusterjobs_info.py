 # -*- coding: utf-8 -*-
"""
Spyder Editor

This is a script file to summarise info of cluster jobs
"""



import os
import glob
import pandas as pd

# on the computer farm - yalan/NNRtest/CAD_UKBIOBANK_200921152632
starttime = '200921152632'

def get_info_from_stdout(cluster, starttime):
	"""
	get some information from stdout files, including: 'run_time', 'cpu_time', 'max_mem', 'ave_mem', and 'n_genes'
	"""
	f = open('CAD_UKBIOBANK-' + cluster + '_' + starttime + '.out')
	outtext = f.readlines()
	f.close()
	#failed_fm = any('failed finemapping\n' in x for x in outtext)
	#n_genes = float("NaN") if failed_fm else [x for x in outtext if 'genes associated\n' in x][0].split(' ')[-3]
	#n_genes = float("NaN") if [x for x in outtext if 'failed finemapping\n' in x] != [] else [x for x in outtext if 'genes associated\n' in x][0].split(' ')[-3]
	return [cluster, outtext[-11].split()[-2], outtext[-19].split()[-2], outtext[-18].split()[-2] if 'MB' in outtext[-18] else 'N/A', outtext[-17].split()[-2] if 'MB' in outtext[-17] else 'N/A'] #, failed_fm, n_genes

sumdf = pd.concat([pd.DataFrame([get_info_from_stdout(cluster, starttime)]) for cluster_folder in glob.glob('CAD_UKBIOBANK-chr*_clusters') for cluster in os.listdir(cluster_folder)], ignore_index=True)
sumdf.columns = ['cluster', 'run_time', 'cpu_time', 'max_mem', 'ave_mem']#, 'failed_fm', 'n_genes'
sumdf.to_csv('CAD_UKBIOBANK_ks1_km2_spcl_jobsum_' + starttime + '.txt', header=True, index=False, sep='\t', na_rep='N/A')
