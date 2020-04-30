# Analysing your own summary statistics with parallelisation

## Pre-checking

In parallelpostgap.py, it will check several important features of GWAS summary statistics uploaded by users before carrying out the analysis.

1. What is the type of the file?
	* Excel files will be converted into tab-separated files;
	* If the suffix of filename is ".out", it will be changed to ".tsv" to avoid conflictting with the output file.
2. Is the file tab-delimited?
	* Tab-delimited is optimal;
	* If the file is separated by some common sparators, including space, comma, semicolon, colom and pipe, it can be changed to tab-separated automatically.
3. It must have the following columns with assigned column names:
   * variant_id
   * p-value
   * beta

If the GWAS summary statistics fails to meet the last two requirements, it will give errors. Please make sure your dataset is tab-delimited (optimal) and has correct column names beforehand.

## Running

Please specify several parameters, when running parallelpostgap.py.

* --database_dir (optional): the directory where data files are stored, and the default directory is "databases/";
* --summary_stats (required): the location of your GWAS summary statistics;
* --time_interval (optional): the sleeping time (in second) added between analyses of chromosomes, which is 180s by default.

The command to run it:

`python scripts/precheck/parallelpostgap.py --summary_stats tests/sample_data/example_GWASsumstats.xlsx --time_interval 30`

GWAS dataset will be splitted into small datasets by chromosomes to run analyses in parellel.

Here we use one GWAs summary statistics to show how to run the command. *Summary statistics were downloaded from the NHGRI-EBI GWAS Catalog (Buniello, MacArthur et al., 2019) for study GCST009246 (Liu C et al., 2019) downloaded on 27/04/2020.* One column name was changed from "SNP" to "variant_id" to carry out analysis.

## Outputs

Currently, Bayesian mode is on by default, and two outputs will be produced and saved under the same directory as summary statistics.

1. One file with the suffix "\_results.txt" is the default output.
2. The other file with the suffix "\_output2.txt" is the output generated in Bayesian mode.
