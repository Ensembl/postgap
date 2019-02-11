POSTGAP
=======

Copyright holder: [EMBL-European Bioinformatics Institute](http://www.ebi.ac.uk) (Apache 2 License)

This script is designed to automatically finemap and highlight the causal variants behind GWAS results by cross-examining GWAS, population genetic, epigenetic and cis-regulatory datasets.

Its original design was based on [STOPGAP](http://www.nature.com/ng/journal/v47/n8/full/ng.3314.html). It takes as input a disease identifier, extracts associated SNPs via GWAS databases, expands them by LD, then searches an array of regulatory and cis-regulatory databases for gene associations.

![Pipeline diagram](https://github.com/Ensembl/postgap/blob/master/POSTGAP%20pipeline.png "Pipeline diagram")

Installing
----------

Add the ```lib/``` directory to your ```$PYTHONPATH``` environment variable.

Installing dependencies
-----------------------

To install all dependencies run ```sh install_dependencies.sh```

Add the ```bin``` directory to your ```$PATH``` environment variable.

Flatfile preparation
--------------------

* Via the FTP site (*recommended*)

  The following script downloads a bunch of files into PWD.
  ```sh download.sh```

  Ideally, save these files in a separate directory, which we will call ```databases_dir````

  Everytime you run POSTGAP, add ```--database_dir /path/to/databases_dir``` to the command line.

* Manually (*sloooow*)
  This script will create a ```databases_dir``` directory for you:
  1. Type ```make download``` to download public databases.
  2. Type ```make process``` to preprocess the databases. **Warning** this may take days as it needs to split the entire 1000 Genomes files by population.

Running
-------

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

Direct data upload
------------------

To short cut the GWAS databases and enter you own data with a file:
```
python POSTGAP.py --summary_stats my_stats.txt
```

The summary statistics file should be tab delimited with the following columns:
* Chromosome (GRCh37)
* Position (GRCh37)
* MarkerName
* Effect_allele
* Non_Effect_allele
* Beta
* SE
* Pvalue

Bayesian mode (EXPERIMENTAL)
----------------------------

For an EFO, you can trigger the Bayesian calculations with:

```
python POSTGAP.py --efos EFO_0000196 --bayesian
```

In this case, POSTGAP produces an output file, 'postgap_output', which can be displayed as:
```
python extract_data.py
```

Output
------

By default, the script writes out a tab delimited file to standard out.

If you wish, you can redirect this into a file:

```
python POSTGAP.py --disease autism --output results.txt
```

or an SQLITE database:

```
python POSTGAP.py --disease autism --db results.sqlite3
```

If you want a JSON dump of all the data retrieved by the pipeline:

```
python POSTGAP.py --disease autism --output results.json --json
python POSTGAP.py --disease autism --json
```

Testing
-------

```
python postgap_and_tests.py --database_dir /path/to/postgap/databases --efos EFO_0008263 --output EFO_0008263.txt --GWAS GWAS_Catalog
```

More Info
---------

Checkout out our [Wiki](https://github.com/Ensembl/postgap/wiki)
