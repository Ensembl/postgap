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

  Everytime you run POSTGAP, add ```--databases_dir /path/to/databases_dir``` to the command line.

* Manually (*sloooow*)
  This script will create a ```databases_dir``` directory for you:
  1. Type ```make download``` to download public databases.
  2. Type ```make process``` to preprocess the databases. **Warning** this may take days as it needs to split the entire 1000 Genomes files by population.

Running
-------

By default, run from the root directory the command: 

```
python scripts/POSTGAP.py --disease autism  
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

Output
------

By default, the script writes out a tab delimited file to standard out.

If you wish, you can redirect this into a file:

```
python scripts/POSTGAP.py --disease autism --output resutls.txt
```

or an SQLITE database:

```
python scripts/POSTGAP.py --disease autism --db results.sqlite3
```

If you want a JSON dump of all the data retrieved by the pipeline:

```
python scripts/POSTGAP.py --disease autism --output resutls.txt --json
python scripts/POSTGAP.py --disease autism --json
```

Testing
-------

```
cat scripts/testing/all_efos.txt | xargs -n1 python scripts/POSTGAP.py --efos > table.tsv
```

