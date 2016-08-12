CTTV 24 Project
===============

Copyright holder: [EMBL-European Bioinformatics Institute](http://www.ebi.ac.uk) (Apache 2 License)

This script is designed to automatically finemap and highlight the causal variants behind GWAS results by cross-examining GWAS, population genetic, epigenetic and cis-regulatory datasets.

Its original design was based on [STOPGAP](). 

Installation
------------

Type '''make''' to download public databases.

This will create a "databases" directory which you can move anywhere. Save the location of that directory at the top of the gwas_to_genes.pl script. 

Running
-------

Must specify 1000 Genomes population to calcuate LD from. Must be indexed bcfs
By default, run with the command: 

'''
python gwas_to_genes.py --disease autism --database_dir ~/lustre2/CTTV24/databases --populations_dir /hps/nobackup/stegle/users/horta/dataset/1000G/LD_calculator/data/processed

'''

Multiple disease names can be provided.

Optional parameters:
--databases\ $database_dir : manually set the location of the database directory
