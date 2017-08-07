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

import argparse
import pickle

def main():
	'''
	
		Loads a pickle file, stringifies and prints it.
		
		Useful when investigating pickled objects that create summaries of themselves like the finemap objects.
		
		E.g.:
		
			python scripts/stringify_pickle_object.py --pickle_file /hps/nobackup/production/ensembl/mnuhn/postgap/work_dir_for_diabetes8/diabetes/gwas_clusters/gwas_cluster_with_247_snps_around_rs1727313/cis_regulatory_evidence_from_eqtl_201_snps_linked_to_ENSG00000111328_in_Whole_Blood.pickle
	
	'''
	
	parser = argparse.ArgumentParser(description='Load and print')
	parser.add_argument('--pickle_file')
	
	options = parser.parse_args()

	object_list = []

	with open(options.pickle_file, 'rb') as pickle_fh:
		try:
			while True:
				current_object = pickle.load(pickle_fh)
				object_list.append(current_object)
		except EOFError:
			print "Loaded " + str(len(object_list)) + " objects from " + options.pickle_file

	from postgap.Summarisers import summarise
	
	for obj in object_list:
		print(summarise(obj));

if __name__ == "__main__":
	main()

