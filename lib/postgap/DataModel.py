#! /usr/bin/env python

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

import collections

Region = collections.namedtuple("Region", ['chrom', 'start', 'end'])
SNP = collections.namedtuple(
	"SNP", 
	[
		'rsID',
		'chrom',
		'pos',
		'approximated_zscore'
	]
)
Gene = collections.namedtuple("Gene", ['name', 'id', 'chrom', 'tss', 'biotype'])
Disease = collections.namedtuple('Disease', ['name', 'efo'])
GWAS_Association = collections.namedtuple(
	'GWAS_Association', 
	[
		'snp',
		'disease',
		'reported_trait',
		'pvalue',
		'pvalue_description',
		'sample_size',
		'source',
		'publication',
		'study',
		'odds_ratio',
		'odds_ratio_ci_start',
		'odds_ratio_ci_end',
		'beta_coefficient',
		'beta_coefficient_unit',
		'beta_coefficient_direction',
		# rest_hash is no longer needed, risk_alleles_present_in_reference 
		# gets set when the object is created and that is all that is ever 
		# needed.
		'rest_hash',
		'risk_alleles_present_in_reference'
	]
)
GWAS_SNP = collections.namedtuple(
	'GWAS_SNP', 
	[
		'snp',
		'pvalue',
		'z_score',
		'evidence',
		'beta'
	]
)

GWAS_Cluster = collections.namedtuple(
	'GWAS_Cluster', 
	[
		'gwas_snps',
		'ld_snps',
		'ld_matrix',
		'z_scores',
		'betas',
		'mafs',
		'annotations',
		'gwas_configuration_posteriors',
		'lambdas'
	]
)

Cisregulatory_Evidence = collections.namedtuple(
	'Cisregulatory_Evidence', 
	[
		'snp',
		'gene',
		'score',
		'source',
		'study',
		'tissue',
		'info',
		'z_score',
		'pvalue',
		'beta'
	]
)
Regulatory_Evidence = collections.namedtuple('Regulatory_Evidence', ['snp','score','source','study','tissue','info'])
GeneSNP_Association = collections.namedtuple('GeneSNP_Association', ['gene', 'snp', 'score', 'rank', 'intermediary_scores', 'cisregulatory_evidence', 'regulatory_evidence'])
GeneCluster_Association = collections.namedtuple('GeneCluster_Association', ['gene', 'cluster', 'score', 'collocation_posterior', 'evidence', 'eQTL_hash', 'r2'])
