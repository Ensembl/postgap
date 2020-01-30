# Post-GWAS Analysis Pipeline

Copyright holder: [EMBL-European Bioinformatics Institute](http://www.ebi.ac.uk) (Apache 2 License)

This script is designed to automatically finemap and highlight the causal variants behind GWAS results by cross-examining GWAS, population genetic, epigenetic and cis-regulatory datasets.

Its original design was based on [STOPGAP](http://www.nature.com/ng/journal/v47/n8/full/ng.3314.html). It takes as input a disease identifier, extracts associated SNPs via GWAS databases, expands them by LD, then searches an array of regulatory and cis-regulatory databases for gene associations.

![Pipeline diagram](https://github.com/Ensembl/postgap/blob/master/POSTGAP%20pipeline.png "Pipeline diagram")

# Installation

## Virtual machine

If you wish to shortcut all of the instructions below, you can simply use our [VirtualBox](https://www.virtualbox.org/)  [virtual machine](http://ftp.ebi.ac.uk/pub/databases/post-gwas/postgap.ova).

## Installing the Python library

Add the ```lib/``` directory to your ```$PYTHONPATH``` environment variable.

## Installing dependencies

The `scripts/installation/ubuntu_environment.sh` describes a recipe to install all basic C and Python dependencies on a fresh ubuntu server (requires root access). 

To install all binformatic dependencies run ```sh scripts/installations/install_dependencies.sh```. 

Add the ```./bin/``` directory to your ```$PATH``` environment variable.

## Flatfile preparation

### Via the FTP site (*recommended*)

  The following script downloads a bunch of files into ```$PWD```:
  ```
  sh scripts/installation/download.sh
  ```

  Ideally, save these files in a separate directory, which we will call ```databases_dir```.

### Manually (*sloooow*)
  The following will create a ```databases_dir``` directory for you:
  ```
  cd scripts/build_data_files
  make download
  make process
  ```
  **Warning** this may take days as it needs to split the entire 1000 Genomes files by population.

# Running
 
Every time you run POSTGAP, add ```--database_dir /path/to/databases_dir``` to the command line, where the database directory path corresponds to the directory created above.

By default, run from the root directory the command: 

```
python POSTGAP.py --disease autism  
```

Multiple disease names can be provided.

You can also provide a list of EFOs:

```
python POSTGAP.py --efos EFO_0000196
```

Or an rsID:

```
python POSTGAP.py --rsID rs10009124
```

Or a manually defined variant:

```
python POSTGAP.py --coords my_variant 1 1234567 
```

## Analysing your own summary statistics

To short cut the GWAS databases and enter you own data with a file:
```
python POSTGAP.py --summary_stats tests/sample_data/example.tsv
```

The summary statistics file should be tab delimited and follow the [GWAS Catalog recommentations](https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics).

In particular, it must have the following columns:
- variant_id
- p-value
- beta

## Bayesian mode (EXPERIMENTAL)

For an EFO, you can trigger the Bayesian calculations with:

```
python POSTGAP.py --efos EFO_0000196 --bayesian --output2 output2.txt
```

In this case, POSTGAP produces a tab-delimited output file, 'output2.txt'. The columns represent:
1. Gene ID
2. Cluster description
3. SNP ID
4. Colocalisation posterior probability at that SNP
5. Tissue
6. Colocalisation posterior probability over the whole cluster

It can be displayed as:
```
python scripts/present_results/postgap_html_report.py --result_file output2.txt --template scripts/present_results/geneReport.html --output report.html
```

### Output

By default, the script writes out a tab delimited file to standard out.

If you wish, you can redirect this into a file:

```
python POSTGAP.py --disease autism --output results.txt
```

If you want a JSON dump of all the data retrieved by the pipeline:

```
python POSTGAP.py --disease autism --output results.json --json
python POSTGAP.py --disease autism --json
```

# Testing

You can check the output with the following commands using the [data tests](./tests/README.md).

# More Info

Check out our [Wiki](https://github.com/Ensembl/postgap/wiki)
