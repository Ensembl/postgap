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
from postgap.DataModel import Region

DATABASES_DIR = None
SPECIES = None
DEBUG = True

BLACKLISTED_REGIONS = [
	Region(chrom = '6', start=28477797, end=33448354) # MHC
]

EVIDENCE_WEIGHTS = {
	'Regulome': 1,
	'VEP': 1,
	'GTEx': 1,
	'Fantom5': 1,
	'DHS': 1,
	'PCHIC': 1,
	'nearest_gene': 1
}

# These get set in POSTGAP.py
work_directory = "notset"

finemap_gwas_clusters_directory = None
finemap_eqtl_clusters_directory = None
