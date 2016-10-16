CTTV 24 Project
===============

Copyright holder: [EMBL-European Bioinformatics Institute](http://www.ebi.ac.uk) (Apache 2 License)

This script is designed to automatically finemap and highlight the causal variants behind GWAS results by cross-examining GWAS, population genetic, epigenetic and cis-regulatory datasets.

Its original design was based on [STOPGAP](). 

Installation
------------

Type '''make download''' to download public databases.
Type '''make process_databases''' to preprocess the databases. 
Type ''' make process_1000G to preprocess 1000Genomes data from the raw vcf.gz files. This can take a long time, it is recommended to download the processed files directly.

Running
-------

By default, run from the root directory the command: 

'''
python scripts/gwas_to_genes.py --disease autism  
'''

Multiple disease names can be provided.

