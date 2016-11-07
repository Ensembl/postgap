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
def concatenate(list):
	"""

		Shorthand to concatenate a list of lists
		Args: [[]]
		Returntype: []

	"""
	return sum(filter(lambda elem: elem is not None, list), [])

def concatenate_hashes(list):
	"""

		Shorthand to concatenate a list of lists
		Args: [[]]
		Returntype: []

	"""
	return dict(sum(map(lambda X: X.items(), filter(lambda elem: elem is not None, list)), []))

def chunks(l, n):
	for i in range(0, len(l), n):
		yield l[i:i+n]

