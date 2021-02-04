# -*- coding: utf-8 -*-
'''
Spyder Editor

This is a script file to generate heatmaps, to visualise genes associated with one cluster in tissues.
a plot for each cluster. axes are genes and tissues, and the color represents cluster_CLPP.
'''



import os
import pandas as pd
from scipy import optimize
import matplotlib.pyplot as plt

starttime = '210201173756'

os.chdir('C:/Users/yalan/Documents/POSTGAP/forNovonordisk/test_' + starttime + '/')

# load in output2 file
df_op2 = pd.read_csv('CAD_UKBIOBANK_output2_ks1_km2_spcl_' + starttime + '.txt', delimiter='\t', header=0, low_memory=False)
df_op2.columns = ['gene_id', 'cluster', 'rsID', 'CLPP_SNP', 'tissue', 'CLPP_cluster']

# get a list of unique clusters
cluster_list = df_op2['cluster'].unique()

# summarise info: clusters, n_SNPs, max_cluster_CLPPs
df = pd.DataFrame([[c, len(df_op2.loc[df_op2['cluster'] == c, 'rsID'].unique()), df_op2.loc[df_op2['cluster'] == c, 'CLPP_cluster'].max()] for c in cluster_list], columns=['cluster', 'size', 'max_CLPP'])
df = df.sort_values('size')

# fit trend curve
def test_func(x, a):
	return a / x

params, params_covariance = optimize.curve_fit(test_func, df['size'], df['max_CLPP'])

# make a plot to show the correlation between size of clusters and max cluster_CLPPs
fig = plt.figure()
plt.scatter(df['size'], df['max_CLPP'], s=16, label='Data')
plt.xlabel('Number of SNPs in the cluster')
plt.ylabel('Maximum CLPP_cluster')
#plt.plot(df['size'], test_func(df['size'], params[0]), color='tab:cyan', label='Fitted function: ' + str(round(params[0], 2)) + ' / x')
plt.plot(df['size'], test_func(df['size'], 1), color='tab:cyan', label='Fitted function: 1 / x')
plt.title('max_iter = 100 & kstart = 1 & kmax = 2')
plt.legend(loc='best')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
fig.set_size_inches(10, 6)
fig.savefig('clusters10+_size_vs_maxCLPP_a1.png', dpi=600, transparent=False, bbox_inches='tight')
plt.close()
