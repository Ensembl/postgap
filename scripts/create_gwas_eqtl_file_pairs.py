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
import logging
import logging.config

logging.config.fileConfig('configuration/logging.conf')
logger = logging.getLogger(__name__)

def main():
	'''
	
python scripts/create_gwas_eqtl_file_pairs.py \
    --gwas_directory              finemap/EFO_0000203/gwas_clusters \
    --eqtl_directory              finemap/EFO_0000203/eqtl_clusters \
    --output_posteriors_directory finemap/EFO_0000203/posteriors \
    --output_commands_file        finemap/EFO_0000203/compute_all_posteriors.bash

	'''
	
	import postgap.Globals
	postgap.Globals.DATABASES_DIR = '/nfs/nobackup/ensembl/mnuhn/postgap/databases/'

	parser = argparse.ArgumentParser(description='Run finemap')
	parser.add_argument('--gwas_directory')
	parser.add_argument('--eqtl_directory')
	parser.add_argument('--output_posteriors_directory')
	parser.add_argument('--output_commands_file')
	
	options = parser.parse_args()

	gwas_directory              = options.gwas_directory
	eqtl_directory              = options.eqtl_directory
	output_posteriors_directory = options.output_posteriors_directory
	output_commands_file        = options.output_commands_file

	logger.info( "gwas_directory              = " + gwas_directory )
	logger.info( "eqtl_directory              = " + eqtl_directory )
	logger.info( "output_posteriors_directory = " + output_posteriors_directory )
	logger.info( "output_commands_file        = " + output_commands_file )
	
	gwas_files = ls_files_from_directory(gwas_directory)
	eqtl_files = ls_files_from_directory(eqtl_directory)
	
	import json
	logger.info( "Got gwas cluster files: " + json.dumps(gwas_files) )
	logger.info( "Got eqtl files: "         + json.dumps(eqtl_files) )
	
	finemap_run_parameters = collections.namedtuple(
		'finemap_run_parameters', [
			'gwas_cluster_file', 
			'eqtl_cluster_file',
			'posteriors_file'
		]
	)
	
	finemap_run_parameters_list = []
	
	for gwas_file in gwas_files:
		for eqtl_file in eqtl_files:
			
			import os
			gwas_file_base = os.path.splitext(gwas_file)[0]
			eqtl_file_base = os.path.splitext(eqtl_file)[0]
			
			posteriors_file = gwas_file_base + "__" + eqtl_file_base + ".pickle"
			
			current_finemap_run_parameters = finemap_run_parameters(
				gwas_cluster_file = gwas_directory              + "/" + gwas_file,
				eqtl_cluster_file = eqtl_directory              + "/" + eqtl_file,
				posteriors_file   = output_posteriors_directory + "/" + posteriors_file
			)
			
			#from pprint import pformat
			#logger.info( "Parameter set: " + pformat(current_finemap_run_parameters))
			
			finemap_run_parameters_list.append(current_finemap_run_parameters)
	
	run_finemap_cmd = "python scripts/run_finemap.py --gwas_cluster_file %s --eqtl_cluster_file %s --output_posteriors_file %s"
	
	output_commands_fh = open(output_commands_file, 'w')
	
	for finemap_run_parameters in finemap_run_parameters_list:
		
		command = run_finemap_cmd % (
			finemap_run_parameters.gwas_cluster_file,
			finemap_run_parameters.eqtl_cluster_file,
			finemap_run_parameters.posteriors_file
		
		)
		logger.info("Generated: " + command)
		output_commands_fh.write(command)
		output_commands_fh.write("\n")
	
	output_commands_fh.close()
	logger.info("Commands have been written to: " + output_commands_file)
	
	logging.info("All done.");

def ls_files_from_directory(directory):
	
	from os import listdir
	from os.path import isfile, join
	
	onlyfiles = [f for f in listdir(directory) if isfile(join(directory, f))]
	onlyfiles.sort()
	
	return onlyfiles

if __name__ == "__main__":
	main()

