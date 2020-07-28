# Analysing your own summary statistics with parallelisation

## Pre-checking

In parallelpostgap.py, it will check several important features of GWAS summary statistics uploaded by users before carrying out the analysis.

1. What is the type of the file?
	* Excel files will be converted into tab-separated files;
	* If the suffix of filename is ".out", it will be changed to ".tsv" to avoid conflictting with the output file.
	* If no suffix in the filename, it will use '.tsv'.
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
* --**summary_stats** (required): the location of your GWAS summary statistics;
* --split_cluster (optional): whether to analyze in parallel by splitting by clusters, and False by default;
* --submit_interval (optional): the sleeping time (in second) added before starting the next analyses of chromosomes / clusters, which is 90s by default;
* --check_interval (optional): the sleeping time (in second) added between checks whether all the analyses are done, which is 300s by default.
parser.add_argument('--check_interval', type=int, default=300, help='sleeping time (s) added to check again if all the analyses of chromosomes have been done')
* --memory (optional): the memory of a node (in GB) requested to carry out analyses, which is 2GB by default;
* --kstart (optional): the minimum number of possile causal variants in each genetic cluster assumed, which is 1 by default;
* --kmax (optional): the maximum number of possile causal variants in each genetic cluster assumed, which is 1 by default.

The command to run it:

`python scripts/precheck/parallelpostgap.py --summary_stats tests/sample_data/example_GWASsumstats.xlsx --submit_interval 15 --check_interval 900 --memory 3 --kstart 1 --kmax 2 --split_cluster`

Here we use one GWAs summary statistics to show how to run the command. *Summary statistics were downloaded from the NHGRI-EBI GWAS Catalog (Buniello, MacArthur et al., 2019) for study GCST009246 (Liu C et al., 2019) downloaded on 27/04/2020.* One column name was changed from "SNP" to "variant_id" to carry out analysis.

With the parameter setting above, it will separate the input GWAS summary statistics both by chromosomes and by clusters, and carry out analyses in parallel. It will request an initial memory of 3GB for analyses, which can be increased up to 5x in case it's not big enough. Analyses will be started with a time gap of 15s, and will be checked every 900s until all are finished. And it's assumed that there would be at least 1 causal variant and at most 2 causal variants at each genetic cluster.


## Outputs

Currently, Bayesian mode is on by default, which means finemapping will be performed. Two outputs will be produced and saved under the same directory as summary statistics.

1. One file with the suffix "\_results.txt" is the default output.
2. The other file with the suffix "\_output2.txt" is the output generated in Bayesian mode.
