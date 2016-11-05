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

import Globals
from DataModel import *
import BedTools

class Reg_source(object):
	def run(self, ld_snps, tissues):
		"""

			Extract score at sns of interest
			Args:
			* [ SNP ]
			* [ string ]
			Returntype: [ Regulatory_Evidence ]

		"""
		assert False, "This stub should be defined"

class Regulome(Reg_source):
	display_name = "Regulome"
	def run(self, ld_snps, tissues):
		"""

			Extract Regulome score at sns of interest
			Args:
			* [ SNP ]
			* [ string ]
			Returntype: [ Regulatory_Evidence ]

		"""
		snp_hash = dict( (snp.rsID, snp) for snp in ld_snps)
		intersection = BedTools.overlap_snps_to_bed(ld_snps, Globals.DATABASES_DIR + "/Regulome.bed")
		res = filter (lambda X: X.score, (self.get_regulome_evidence(feature, snp_hash) for feature in intersection))

		if Globals.DEBUG:
			print "\tFound %i regulatory variants in Regulome" % (len(res))

		return res

	def get_regulome_evidence(self, feature, snp_hash):
		"""

			Extract Regulome score from bedtools output
			Args:
			* string
			Returntype: Regulatory_Evidence

		"""

		'''
		First 4 columns are from Regulome file, the next 4 are LD SNP coords
		Regulome file format:
		1. chrom
		2. start
		3. end
		4. category

		LD SNP coords:
		5. chrom
		6. start
		7. end
		8. rsID
		'''
		return Regulatory_Evidence(
				snp = snp_hash[feature[7]],
				source = self.display_name,
				score = 2 if feature[3][0] == '1' or feature[3][0] == '2' else 1,
				study = None,
				tissue = None
			)

sources = Reg_source.__subclasses__()
