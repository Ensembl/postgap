# Analysing your own summary statistics

## Pre-checking

In prepareGWASsumstat.py, it will check two important features before carrying out postgap.py.

1. The summary statistics file should be tab delimited;
2. it must have the following columns with assigned column names:
   * variant_id
   * p-value
   * beta

If the GWAS summary statistics doesn't meet these two requirements, it will give errors. Please make sure your dataset is tab delimited and has correct column names beforehand.

## Running

Please specify the filename of your GWAS summary statistics when running prepareGWASsumstat.py.

`python scripts/precheck/prepareGWASsumstat.py --filename Shrine.txt`

GWAS dataset will be splitted into small datasets by chromosomes to run analysis in parellel.

## Outputs

Currently, Bayesian mode is on by default, and two outputs will be produced. One file with the suffix "\_res.txt" is the default output. The other file with the suffix "\_op2.txt" is the output generated in Bayesian mode. The columns represent:

1. Gene ID
2. Cluster description
3. SNP ID
4. Colocalisation posterior probability at that SNP
5. Tissue
6. Colocalisation posterior probability over the whole cluster
