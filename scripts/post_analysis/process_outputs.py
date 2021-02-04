# -*- coding: utf-8 -*-
'''
Spyder Editor

This is a script file to process the outputs of an analysis to create the dataframe for making plots
'''



import os
import pandas as pd

# get command line arguments
parser = argparse.ArgumentParser(description='to read command line arguments')
parser.add_argument('--starttime', type=str, default='210130231014', help='the time when the analysis was carried out')

args = parser.parse_args()
starttime = args.starttime

os.chdir('test_' + starttime + '/')
#os.chdir('C:/Users/yalan/Documents/POSTGAP/forNovonordisk/test_' + starttime + '/')

# load in results & output2 files
df_res = pd.read_csv('CAD_UKBIOBANK_results_ks1_km2_spcl_' + starttime + '.txt', delimiter='\t', header=0, low_memory=False)
df_op2 = pd.read_csv('CAD_UKBIOBANK_output2_ks1_km2_spcl_' + starttime + '.txt', delimiter='\t', header=0, low_memory=False)
df_op2.columns = ['gene_id', 'cluster', 'rsID', 'CLPP_SNP', 'tissue', 'CLPP_cluster']

# list SNP - associated gene pairs in results
snp_gene_pairs = df_res[['ld_snp_rsID', 'gene_id']].drop_duplicates()

# refine output2: remove SNP - gene pairs that in output2 but not in results
df_plot = df_op2.loc[[any((snp_gene_pairs['ld_snp_rsID'] == row['rsID']) & (snp_gene_pairs['gene_id'] == row['gene_id'])) for index, row in df_op2.iterrows()]]
df_plot.to_csv('CAD_UKBIOBANK_df_plots_' + starttime + '.tsv', sep='\t', header=True, index=False)

# run under the directory yalan/NNRtest/
#gsub -m 50 -n df_plot -d . -q production-rh74 'python process_outputs.py --starttime 210130231014' -r
