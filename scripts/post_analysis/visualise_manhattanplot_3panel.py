# -*- coding: utf-8 -*-
'''
Spyder Editor

This is a script file to generate gene-tissue specific 3-panel Manhattan plots.
	a 3-panel plot for all the clusters associated with a gene under a tissue;
	axes are position of SNPs and SNP-specific GWAS p-value, eQTL p-value and SNP CLPP.
'''



import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

starttime = '201130093846'
path_to_chr_files = 'C:/Users/yalan/Documents/POSTGAP/forNovonordisk/CAD_UKBIOBANK_chr/' #GWAS summary statistics of chromosomes

os.chdir('C:/Users/yalan/Documents/POSTGAP/forNovonordisk/test_' + starttime + '/')

def create_manhattan(suminfo, data, gene_list):
	gene = suminfo['gene_id']
	tissue = suminfo['tissue']
	
	colorlist = ['b', 'g', 'r', 'c', 'm', 'y']
	
	# create the 3-panel Manhattan plot
	fig, (ax1, ax2, as3) = plt.subplots(3, sharex=True, sharey=False)
	fig.suptitle(gene_list[gene] + ' (' + ('' if suminfo['found_in_paper'] else 'not ') + 'found in paper) in ' + tissue)
	ax1.scatter(data['pos'], -np.log10(data['pv_GWAS']), c=[colorlist[i % len(colorlist)] for cl in data['cluster'] for i, c in enumerate(data['cluster'].unique()) if cl == c], s=8)
	ax1.set_ylabel('-log10(pv_GWAS)', fontsize=8)
	ax2.scatter(data['pos'], -np.log10(data['pv_eQTL']), c=[colorlist[i % len(colorlist)] for cl in data['cluster'] for i, c in enumerate(data['cluster'].unique()) if cl == c], s=8)
	ax2.set_ylabel('-log10(pv_eQTL)', fontsize=8)
	as3.scatter(data['pos'], data['CLPP_SNP'], c=[colorlist[i % len(colorlist)] for cl in data['cluster'] for i, c in enumerate(data['cluster'].unique()) if cl == c], s=8)
	as3.set_ylabel('CLPP_SNP', fontsize=8)
	plt.xlabel('chr ' + chr)
	fig.set_size_inches(8, 8)
	fig.savefig('ManhattanPlots/' + gene_list[gene] + '-' + tissue + '_MP.png', dpi=600, transparent=False, bbox_inches='tight')
	plt.close()



# load in results & df_plot files
df_res = pd.read_csv('forNovonordisk/test_' + starttime + '/CAD_UKBIOBANK_results_ks1_km2_spcl_' + starttime + '.txt', delimiter='\t', header=0, low_memory=False)
df_plot = pd.read_csv('forNovonordisk/test_' + starttime + '/CAD_UKBIOBANK_output_plots.tsv', delimiter='\t', header=0, low_memory=False)

# get gene_id - gene_symbol pairs in results, and add to df_plot
gene_symbol_pairs = df_res[['gene_id', 'gene_symbol']].drop_duplicates()
gene_list = dict(zip(gene_symbol_pairs.gene_id, gene_symbol_pairs.gene_symbol))

# update df_plot
df_plot = df_plot.reindex(columns=['gene_id', 'tissue', 'cluster', 'CLPP_cluster', 'rsID', 'chr', 'pos', 'CLPP_SNP', '-log10_CLPP_SNP', 'pv_eQTL', '-log10_pv_eQTL', 'pv_GWAS', '-log10_pv_GWAS'], fill_value=0)

for i, row in df_plot.iterrows():
	refrow = df_res.loc[(df_res['ld_snp_rsID'] == row['rsID']) & (df_res['ld_snp_is_gwas_snp'] == 1) & (df_res['gene_id'] == row['gene_id'])]
	assert len(refrow) == 1, 'debug here!'
	# get LD SNPs' chr and pos, GWAS p-value when LD SNP IS GWAS SNP, and get gene-tissue specific eQTL p-value
	df_plot.at[i, 'pv_eQTL'] = 1 - refrow['GTEx_' + row['tissue']].values[0]
	df_plot.at[i, ['chr', 'pos', 'pv_GWAS']] = refrow[['GRCh38_chrom', 'GRCh38_pos', 'gwas_pvalue']].values[0]

df_plot['-log10_CLPP_SNP'] = - np.log10(df_plot['CLPP_SNP'])
df_plot['-log10_pv_eQTL'] = - np.log10(df_plot['pv_eQTL'])
df_plot['-log10_pv_GWAS'] = - np.log10(df_plot['pv_GWAS'])

# summarise relevant info into a table, containing top 10 genes and genes found by both us and the paper (true positive)
df = pd.DataFrame(columns=['cluster', 'gene_id', 'gene_symbol', 'tissue', 'found_in_paper'])
df[['cluster', 'gene_symbol', 'found_in_paper']] = gene_table.loc[:9, ['cluster', 'gene_symbol', 'found_in_paper']].append(gene_table.loc[gene_table['found_in_paper'] == True, ['cluster', 'gene_symbol', 'found_in_paper']], ignore_index=True)
df['gene_id'] = [i for g in df['gene_symbol'] for i in gene_list if gene_list[i] == g]
for i,row in df.iterrows():
	df.at[i, 'tissue'] = df_plot.loc[df_plot.loc[(df_plot['gene_id'] == row['gene_id']) & (df_plot['cluster'] == row['cluster']), 'CLPP_cluster'].idxmax(), 'tissue']

# make Manhattan plots for top 10 genes and genes found by both us and the paper (true positive)
for i, row in df.iterrows():
	gene = row['gene_id']
	tissue = row['tissue']
	
	# get data of clusters associated with this gene under this tissue
	data = df_plot.loc[(df_plot['tissue'] == tissue) & (df_plot['gene_id'] == gene)]
	assert len(data) > 0, 'debug here!'
	#print(str(len(data['cluster'].unique())))
	
	# achieve LD SNPs' chr and pos, GWAS p-value when LD SNP IS GWAS SNP, and get gene-tissue specific eQTL p-value
	data = data.reindex(columns=['gene_id', 'tissue', 'cluster', 'CLPP_cluster', 'rsID', 'chr', 'pos', 'pv_GWAS', 'pv_eQTL', 'CLPP_SNP'], fill_value=0)
	
	data['chr'] = [c.split(':')[0].split('_')[-1] for c in data['cluster']]
	assert len(data['chr'].unique()) == 1, 'clusters from multiple chr under this gene-tissue pair'
	
	chr = data['chr'].unique()[0]
	chr_sum = pd.read_csv(path_to_chr_files + 'CAD_UKBIOBANK-chr' + chr + '.tsv', sep='\t', header=0, low_memory=False)
	
	for k, drw in data.iterrows():
		rsID = drw['rsID']
		refrow = df_res.loc[(df_res['ld_snp_rsID'] == rsID) & (df_res['ld_snp_is_gwas_snp'] == 1) & (df_res['gene_id'] == gene)]
		
		if len(refrow) == 1:
			data.loc[k, 'pv_GWAS'] = refrow['gwas_pvalue'].values[0]
		elif len(refrow) == 0:
			#get pv_GWAS from the chr file
			data.loc[k, 'pv_GWAS'] = chr_sum.loc[chr_sum['variant_id'] == rsID, 'p-value'].values[0]
			
			refrow = df_res.loc[(df_res['ld_snp_rsID'] == rsID) & (df_res['gene_id'] == gene), ['GRCh38_chrom', 'GRCh38_pos', 'GTEx_' + tissue]].drop_duplicates()
			assert len(refrow) == 1, 'ld_snp_is_gwas_snp == 0, & multiple eQTL_pv'
		else:
			assert False, 'debug here!'
		
		# obtain eQTL p-value from output table
		#data.at[k, 'pv_eQTL'] = 1 - refrow['GTEx_' + drw['tissue']].values[0]
		# or temporarily just get it from ensembl VEP
		import requests, sys
		server = 'https://rest.ensembl.org'
		ext = '/eqtl/id/homo_sapiens/' + gene + '?tissue=' + tissue + ';statistic=p-value;variant_name=' + rsID
		r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
		assert r.ok, r.raise_for_status()
		decoded = r.json()
		data.loc[k, 'pv_eQTL'] = decoded[0]['value'] if len(decoded) > 0 else np.nan

		data.at[k, ['chr', 'pos']] = refrow[['GRCh38_chrom', 'GRCh38_pos']].values[0]
	
	# make Manhattan plot
	create_manhattan(row, data, gene_list)
