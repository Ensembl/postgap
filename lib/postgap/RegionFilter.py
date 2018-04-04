
"""

Copyright [1999-2018] EMBL-European Bioinformatics Institute

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
from postgap.DataModel import *
from postgap.Globals import BLACKLISTED_REGIONS

def region_filter(clusters):
	return filter(lambda cluster: not cluster_overlap_regions(cluster, BLACKLISTED_REGIONS), clusters)

def cluster_overlap_regions(cluster, regions):
	return any(cluster_overlap_region(cluster, region) for region in regions)

def cluster_overlap_region(cluster, region):
	return any(snp_overlap_region(snp, region) for snp in cluster.ld_snps)

def snp_overlap_region(snp, region):
	return snp.chrom == region.chrom and snp.pos >= region.start and snp.pos < region.end
