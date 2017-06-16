#!/usr/bin/env python

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

def compute_gwas_clusters_with_finemap_posteriors(gwas_clusters):
	
	from methods.GWAS_Cluster import GWAS_Clusters_ok
	assert(GWAS_Clusters_ok(gwas_clusters))
	
	from methods.GWAS_Cluster import compute_gwas_cluster_with_finemap_posteriors	
	gwas_clusters_with_posteriors = [ 
		compute_gwas_cluster_with_finemap_posteriors(gwas_cluster) 
			for gwas_cluster in gwas_clusters 
	]
	return gwas_clusters_with_posteriors

