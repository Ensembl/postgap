# ------------------------------------------------
# built-ins
import os
import datetime
import itertools

# pipped
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from IPython.display import display, HTML

matplotlib.style.use('ggplot')
# ------------------------------------------------

ID_FIELDS = ['gene_id', 'ld_snp_rsID', 'gwas_snp', 'disease_efo_id', 'gwas_pmid']
ID_FIELD_PAIRS = [
    ['gene_id', 'ld_snp_rsID'],
    ['ld_snp_rsID', 'gwas_snp'],
    ['gwas_snp', 'disease_efo_id'],
    ['disease_efo_id', 'gwas_pmid']
]
G2V_FIELDS = ['VEP', 'Regulome', 'PCHiC', 'GTEx', 'Fantom5', 'DHS', 'Nearest']
V2D_FIELDS = ['gwas_pvalue', 'gwas_beta', 'gwas_odds_ratio', 'gwas_size']
FG_FIELDS = ['PCHiC', 'GTEx', 'Fantom5', 'DHS']
XLIMS_ZERO_ONE_FIELDS = ['Regulome', 'PCHiC', 'GTEx', 'Fantom5', 'DHS', 'Nearest']
STANDARD_FIG_SIZE = (15, 5)
STANDARD_HEAD_NUM = 3

def load_file(filename):
    return pd.read_csv(filename, sep='\t', na_values=['None'])

def print_df(df, index=False):
    display(HTML(df.to_html(index=index)))

def print_hist(series, title='', xlabel='Value', ylabel='Frequency'):
    plt.figure(figsize=STANDARD_FIG_SIZE)
    plt.hist(series, bins=100, log=True)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()

def print_field_hists(pg, field_list, figsize=STANDARD_FIG_SIZE):
    plt.figure(figsize=figsize)
    for (i, c) in enumerate(field_list):
        plt.subplot(len(field_list) // 2 + 1, 2, i + 1)
        plt.hist(pg[c].dropna(), bins=100, log=True)
        plt.title('{} Distribution'.format(c))
        if c in XLIMS_ZERO_ONE_FIELDS:
            plt.xlim(0, 1)
        plt.ylabel('Frequency')
    plt.show()

def calc_g2v_field_hists(pg):
    '''Calculate the distributions of G2V subscores (across unique gene-LD SNP pairs)'''
    # TODO: pd.DataFrame.first() is preferable to pd.DataFrame.nth(0), but 
    #       see issue https://github.com/pandas-dev/pandas/issues/19283.
    #       Upgrade pandas when they have fixed this.
    per_g2v = pg.groupby(['gene_id', 'ld_snp_rsID']).nth(0).reset_index()
    print_field_hists(per_g2v, G2V_FIELDS, figsize=(15, 20))

def calc_v2d_field_hists(pg):
    '''Calculate the distributions of V2D subscores (across unique GWAS SNP-disease pairs)'''
    # TODO: Check this should be per (gwas_snp, disease). Does gwas_study/gwas_pmid etc matter?
    per_v2d = pg.groupby(['gwas_snp', 'disease_efo_id']).nth(0).reset_index()
    print_field_hists(per_v2d, V2D_FIELDS, figsize=(15, 10))

def calc_run_str():
    '''Calculate when and by who the notebook was generated.'''
    now = datetime.datetime.now()
    user = os.environ['USER']
    print('Notebook generated at {} by {}'.format(now.isoformat(), user))

def calc_id_field_counts(pg):
    '''Calculate how many unique values occur for each of ID_FIELDS.'''
    df = pd.DataFrame([[c, pg[c].nunique()] for (i, c) in enumerate(ID_FIELDS)],
                      columns=['field', 'unique_values'])
    print_df(df)

def calc_id_field_max_rows(pg):
    '''Calculate the top few max occurrences for each of ID_FIELDS.'''
    for (i, c) in enumerate(ID_FIELDS):
        groups = pg.groupby(c)
        group_frequencies = groups.size().sort_values(ascending=False)
        df = group_frequencies.to_frame(name='row_occurrences').reset_index(c).head(STANDARD_HEAD_NUM)
        print_df(df)

def calc_field_pair_counts(pg, field_pairs):
    '''Calculate how many unique values occur for each field pair in field_pairs.'''
    df = pd.DataFrame([[c, len(pg.groupby(c).size())] for (i, c) in enumerate(field_pairs)],
                      columns=['field_pair', 'unique_associations'])
    print_df(df)

def calc_g2d_pair_counts(pg):
    calc_field_pair_counts(pg, [['gene_id', 'disease_efo_id']])

def calc_id_field_pair_counts(pg):
    calc_field_pair_counts(pg, ID_FIELD_PAIRS)

def calc_pairwise_degree_dist(pg, field_a, field_b, label_a, label_b):
    '''Calculate the degree distribution across A nodes to B nodes and vice versa.'''
    a_degrees = pg.groupby(field_a)[field_b].nunique()
    b_degrees = pg.groupby(field_b)[field_a].nunique()

    plt.figure(figsize=STANDARD_FIG_SIZE)

    plt.subplot(121)
    plt.hist(a_degrees.values.flatten(), bins=100, log=True)
    plt.title('{}s per {} Distribution'.format(label_b, label_a))
    plt.ylabel('Frequency')
    plt.xlabel('Degree')

    plt.subplot(122)
    plt.hist(b_degrees.values.flatten(), bins=100, log=True)
    plt.title('{}s per {} Distribution'.format(label_a, label_b))
    plt.ylabel('Frequency')
    plt.xlabel('Degree')
    plt.show()

def calc_dist_r2(pg):
    '''Calculate the distribution of r2 (across unique LD SNP-GWAS SNP pairs)'''
    r2 = pg.groupby(['ld_snp_rsID', 'gwas_snp'])['r2'].nth(0)
    print_hist(r2, title='LD (r2) Distribution')

def calc_g2v_field_cross_dists(pg):
    '''Calculate pairwise distributions (as heatmap) of G2V fields.'''
    PAIR = ['gene_id', 'ld_snp_rsID']
    subscore_fields = pg[[*PAIR, *G2V_FIELDS]].groupby(PAIR).nth(0)
    combs = itertools.combinations(G2V_FIELDS, 2)

    plt.figure(figsize=(15, 15))
    for (i, c) in enumerate(combs):
        fx = c[0]
        fy = c[1]
        ix = G2V_FIELDS.index(fx)
        iy = G2V_FIELDS.index(fy)
        l = len(G2V_FIELDS)
        plt.subplot(l, l, (iy * l) + ix + 1)
        plt.hist2d(subscore_fields[fx], subscore_fields[fy], bins=20, norm=LogNorm(), cmap=plt.get_cmap('Reds'))
        if fx in XLIMS_ZERO_ONE_FIELDS:
            plt.xlim(0, 1)
        if fy in XLIMS_ZERO_ONE_FIELDS:
            plt.ylim(0, 1)
        if (ix == 0):
            plt.ylabel('{}'.format(fy))
        if (iy == l - 1):
            plt.xlabel('{}'.format(fx))
    plt.tight_layout()
    plt.show()

def calc_g2v_field_overlap(pg):
    '''Calculate the venn diagram of G2V fields.'''
    per_gene_and_ld_snp = pg.groupby(['gene_id', 'ld_snp_rsID']).nth(0)[G2V_FIELDS]
    g2v_field_venn = per_gene_and_ld_snp.apply(lambda x: str([c for c in x.index if x[c] > 0]), axis=1).value_counts()
    print(g2v_field_venn)
