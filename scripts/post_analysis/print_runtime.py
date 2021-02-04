 # -*- coding: utf-8 -*-
"""
Spyder Editor

This is a script file to quickly know the running time (h) / job status of cluster jobs
"""



import os
import glob

# on the computer farm - yalan/NNRtest/CAD_UKBIOBANK_200921152632
starttime = '200921152632'

# a very quick view of running time (h) / job status
rt = []
for cluster_folder in glob.glob('CAD_UKBIOBANK-chr*_clusters'):
	for cluster in os.listdir(cluster_folder):
		if os.path.exists('CAD_UKBIOBANK-' + cluster + '_' + starttime + '.out'):
			f = open('CAD_UKBIOBANK-' + cluster + '_' + starttime + '.out')
			outtext = f.readlines()
			f.close()

			if 'Successfully completed' in outtext[-23]:
				rt.append(float(outtext[-11].split()[-2])/3600)
			elif 'Exited with' in outtext[-23]:
				print(cluster + 'failed the job')
			else:
				print(cluster + ' unfinished yet')
		else:
			print(cluster + '\'s job stdout has not created')
