#! /usr/bin/env python

"""

Copyright [1999-2016] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License")
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

		 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

"""

"""

	Please email comments or questions to the public Ensembl
	developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

	Questions may also be sent to the Ensembl help desk at
	<http://www.ensembl.org/Help/Contact>.

"""
import sys
import argparse
import collections
import json
import sqlite3
import logging
import logging.config

logging.config.fileConfig('configuration/logging.conf')
logger = logging.getLogger(__name__)

def main():
	'''

python scripts/summarise_finemap_object.py --finemap_pickle_file finemap/EFO_0000203/gwas_posteriors/gwas_cluster_0.pickle

	'''
	
	parser = argparse.ArgumentParser(description='Run finemap')
	parser.add_argument('--finemap_pickle_file')
	
	options = parser.parse_args()

	logger.info( "finemap_pickle_file = " + options.finemap_pickle_file )
	
	import pickle
	pickle_fh       = open(options.finemap_pickle_file, 'rb')
	finemap_object  = pickle.load(pickle_fh)

	logging.info(finemap_object);
	logging.info("All done.");

if __name__ == "__main__":
	main()

