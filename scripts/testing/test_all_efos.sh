#!/bin/bash -Eeu

cat all_efos.txt | xargs -n1 python gwas_to_genes.py --disease
