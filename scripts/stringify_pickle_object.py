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

def main():
	'''
	
		Loads a pickle file, stringifies and prints it.
		
		Useful when investigating pickled objects that create summaries of themselves like the finemap objects.
		
		E.g.:
		
			python scripts/stringify_pickle_object.py --pickle_file finemap/EFO_0000203/joint_posteriors/eqtl_snps_linked_to_ENSG00000168038_in_Adipose_Subcutaneous__gwas_cluster_0.joint_posteriors.pickle
	
	'''
	
	parser = argparse.ArgumentParser(description='Load and print')
	parser.add_argument('--pickle_file')
	
	options = parser.parse_args()

	import pickle
	pickle_fh       = open(options.pickle_file, 'rb')
	finemap_object  = pickle.load(pickle_fh)

	print(str(finemap_object));

if __name__ == "__main__":
	main()

