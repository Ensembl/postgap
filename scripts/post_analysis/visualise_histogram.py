 # -*- coding: utf-8 -*-
"""
Spyder Editor

This is a script file to make histograms of runtime, CPU time and max RAM usage of cluster jobs
"""



import math
#import numpy as np
import os

os.chdir('C:/Users/yalan/Documents/POSTGAP/forNovonordisk/test_' + starttime + '/')

starttime = '201130093846'
sumdf = pd.read_csv('CAD_UKBIOBANK_clusum_ks1_km2_spcl_' + starttime + '.txt', sep='\t', header=0, low_memory=False)


# running time
plt.figure(figsize=(8, 6))
plt.hist(x=sumdf['run_time']/3600, color='olivedrab', alpha=0.7, rwidth=0.85)
#plt.hist(x=np.ceil(sumdf['run_time']/3600), bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
#plt.axvline(x=5, color='r', linestyle='-', linewidth=1)
plt.xlabel('Time(h)')
plt.ylabel('n')
plt.title('Histogram of Running Time of ' + str(sumdf.shape[0]) + ' Clusters')
#plt.text(23, 45, r'$\mu=15, b=3$')
plt.savefig('hist_running_time_' + starttime + '.png', dpi=600, transparent=False, bbox_inches='tight')
plt.close()


# CPU time
plt.figure(figsize=(8, 6))
plt.hist(x=sumdf['cpu_time']/3600, color='olivedrab', alpha=0.7, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
#plt.axvline(x=5, color='r', linestyle='-', linewidth=1)
plt.xlabel('Time(h)')
plt.ylabel('n')
plt.title('Histogram of CPU Time of ' + str(sumdf.shape[0]) + ' Clusters')
#plt.text(23, 45, r'$\mu=15, b=3$')
plt.savefig('hist_cpu_time_' + starttime + '.png', dpi=600, transparent=False, bbox_inches='tight')
plt.close()


# max memory
plt.figure(figsize=(8, 6))
plt.hist(x=sumdf['max_mem']/1000, bins=[i for i in range(math.ceil(sumdf['max_mem'].max()/1000 + 1))], color='olivedrab', alpha=0.7, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
#plt.axvline(x=5, color='r', linestyle='-', linewidth=1)
plt.xlabel('Memory(GB)')
plt.ylabel('n')
plt.title('Histogram of Memory Usage of ' + str(sumdf.shape[0]) + ' Clusters')
#plt.text(23, 45, r'$\mu=15, b=3$')
plt.savefig('hist_max_mem_' + starttime + '.png', dpi=600, transparent=False, bbox_inches='tight')
plt.close()
