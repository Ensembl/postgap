# -*- coding: utf-8 -*-
'''
Spyder Editor

This is a script file to generate heatmaps, to visualise genes associated with one cluster in tissues.
	a plot for each cluster;
	axes are genes and tissues, and the color represents cluster_CLPP.
'''



import os
import pandas as pd
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt

starttime = '201130093846'

os.chdir('C:/Users/yalan/Documents/POSTGAP/forNovonordisk/test_' + starttime + '/')

def create_heatmap(cluster, cluster_data, gene_list):
	cluster_heatmap = cluster_data.pivot(index='gene_id', columns='tissue', values='CLPP_cluster')
	# create the heatmap
	plt.figure()
	sns.set(font_scale=0.5)
	ax = sns.heatmap(cluster_heatmap, linewidths=.3, cmap='YlOrRd', yticklabels=[gene_list[t] for t in cluster_heatmap.index])
	plt.title(cluster, fontsize=12)
	plt.xlabel('Tissues', fontsize=15)
	plt.ylabel('Genes', fontsize=15)
	ax.set_yticklabels(ax.get_yticklabels(), rotation=30)
	fig = ax.get_figure()
	fig.savefig('heatmaps/' + cluster.replace(':', '_') + '_heatmap.png', dpi=200, transparent=False, bbox_inches='tight')
	plt.close()



# load in results and df_plots files
df_res = pd.read_csv('CAD_UKBIOBANK_results_ks1_km2_spcl_' + starttime + '.txt', delimiter='\t', header=0, low_memory=False)
df_plot = pd.read_csv('CAD_UKBIOBANK_df_plot_' + starttime + '.tsv', delimiter='\t', header=0, low_memory=False)

# get gene_id - gene_symbol pairs in results, and add to df_plot
gene_symbol_pairs = df_res[['gene_id', 'gene_symbol']].drop_duplicates()
gene_list = dict(zip(gene_symbol_pairs.gene_id, gene_symbol_pairs.gene_symbol))
#for gene_id in gene_list: df_plot.loc[df_plot['gene_id'] == gene_id, 'gene_symbol'] = gene_list[gene_id]

# get a list of unique clusters from df_plot
cluster_list = df_plot['cluster'].unique()

# generate heatmaps
for cluster in cluster_list:
	cluster_data = df_plot.loc[df_plot['cluster'] == cluster, ['gene_id', 'tissue', 'CLPP_cluster']].drop_duplicates()
	create_heatmap(cluster, cluster_data, gene_list)
