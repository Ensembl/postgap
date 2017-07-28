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
from macpath import basename
from postgap.Integration import cisregulatory_evidence

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

from os.path import * 

def main():
	'''
	
python scripts/create_finemap_pipeline_commands.py \
    --gwas_clusters_directory        /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir_for_diabetes8/diabetes/gwas_clusters/ \
    --output_commands_file /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir_for_diabetes8/diabetes/gwas_clusters/compute_all_posteriors.bash

	'''
	
	import postgap.Globals
	postgap.Globals.DATABASES_DIR = '/nfs/nobackup/ensembl/mnuhn/postgap/databases/'

	parser = argparse.ArgumentParser(description='Run finemap')
	parser.add_argument('--gwas_clusters_directory')
	parser.add_argument('--output_commands_file')
	
	options = parser.parse_args()

	gwas_clusters_directory = options.gwas_clusters_directory
	output_commands_file    = options.output_commands_file

	logger.info( "gwas_clusters_directory = " + gwas_clusters_directory )
	logger.info( "output_commands_file    = " + output_commands_file )
	
	output_commands_fh = open(output_commands_file, 'w')
	
	gwas_files = ls_files_from_directory(gwas_clusters_directory)
		
	for gwas_cluster_file in gwas_files:
		
		(gwas_cluster_file_basename, gwas_cluster_file_extension) = splitext(gwas_cluster_file)
		
		if gwas_cluster_file_extension == ".bash":
			continue
		
		if gwas_cluster_file_extension != ".pickle":
			raise Exception("File " + gwas_cluster_file_extension + " does not end in .pickle")
		
		gwas_cluster_posteriors_output_file = gwas_clusters_directory + "/" + gwas_cluster_file_basename + ".posteriors" + gwas_cluster_file_extension
		
		cmd = create_command_for_gwas_cluster_posteriors(
			gwas_cluster_file                   = gwas_clusters_directory + "/" + gwas_cluster_file, 
			gwas_cluster_posteriors_output_file = gwas_cluster_posteriors_output_file,
		)

		output_commands_fh.write(cmd)
		output_commands_fh.write("\n")
				
		cisregulatory_evidence_directory = gwas_clusters_directory + "/" + gwas_cluster_file_basename
		
		if isdir(cisregulatory_evidence_directory):
			
			cisreg_files = ls_files_from_directory(cisregulatory_evidence_directory)
			
			for cisreg_file in cisreg_files:
				
				(cisreg_file_basename, cisreg_file_extension) = splitext(cisreg_file)
				
				eqtl_posteriors_output_file = cisregulatory_evidence_directory + "/" + cisreg_file_basename + ".posteriors" + cisreg_file_extension
				
				cmd = create_command_for_eqtl_posteriors(
					gwas_cluster_file           = gwas_clusters_directory          + "/" + gwas_cluster_file,
					eqtl_file                   = cisregulatory_evidence_directory + "/" + cisreg_file,
					eqtl_posteriors_output_file = eqtl_posteriors_output_file,  
				)
				
				output_commands_fh.write(cmd)
				output_commands_fh.write("\n")
				
				joint_posteriors_output_file = cisregulatory_evidence_directory + "/" + cisreg_file_basename + "_and_" + gwas_cluster_file_basename + ".joint_posteriors" + cisreg_file_extension
				
				cmd = create_command_for_joint_posteriors(
					gwas_cluster_posteriors_file = gwas_cluster_posteriors_output_file,
					eqtl_posteriors_file         = eqtl_posteriors_output_file,  
					joint_posteriors_file        = joint_posteriors_output_file,
				)
				
				output_commands_fh.write(cmd)
				output_commands_fh.write("\n")
				
	logger.info("Commands have been written to: " + output_commands_file)
	logging.info("All done.");
				

def create_command_for_gwas_cluster_posteriors(gwas_cluster_file, gwas_cluster_posteriors_output_file):
	
	cmd = "python scripts/compute_gwas_posteriors_with_finemap.py \\\n" \
		+ "    --gwas_cluster_file      " + gwas_cluster_file + " \\\n" \
		+ "    --output_posteriors_file " + gwas_cluster_posteriors_output_file + "\n"
	
	return cmd

def create_command_for_eqtl_posteriors(gwas_cluster_file, eqtl_file, eqtl_posteriors_output_file):
	
	cmd = "python scripts/compute_eqtl_posteriors_with_finemap.py \\\n" \
		+ "    --gwas_cluster_file       " + gwas_cluster_file + " \\\n" \
		+ "    --eqtl_cluster_file       " + eqtl_file + " \\\n" \
		+ "    --output_posteriors_file  " + eqtl_posteriors_output_file + "\n"

	return cmd

def create_command_for_joint_posteriors(gwas_cluster_posteriors_file, eqtl_posteriors_file, joint_posteriors_file):
	
	cmd = "python scripts/compute_joint_posteriors.py" \
		+ " --gwas_posterior " 		    + gwas_cluster_posteriors_file + " \\\n" \
		+ " --eqtl_posterior "          + eqtl_posteriors_file + " \\\n" \
		+ " --output_joint_posteriors " + joint_posteriors_file + "\n"

	return cmd

def ls_files_from_directory(directory):
	
	from os import listdir
	from os.path import isfile, join

	
	onlyfiles = [f for f in listdir(directory) if isfile(join(directory, f))]
	onlyfiles.sort()
	
	return onlyfiles

if __name__ == "__main__":
	main()

