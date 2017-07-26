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

import collections

Region = collections.namedtuple("Region", ['chrom', 'start', 'end'])
SNP = collections.namedtuple("SNP", ['rsID', 'chrom', 'pos'])
Gene = collections.namedtuple("Gene", ['name', 'id', 'chrom', 'tss', 'biotype'])
Disease = collections.namedtuple('Disease', ['name', 'efo'])
GWAS_Association = collections.namedtuple('GWAS_Association', ['snp','disease','reported_trait','pvalue','sample_size','source','study'])
GWAS_SNP = collections.namedtuple('GWAS_SNP', ['snp','pvalue', 'evidence'])
GWAS_Cluster = collections.namedtuple('GWAS_Cluster', ['gwas_snps','ld_snps'])
Cisregulatory_Evidence = collections.namedtuple('Cisregulatory_Evidence', ['snp','gene','score','source','study','tissue','info'])
Regulatory_Evidence = collections.namedtuple('Regulatory_Evidence', ['snp','score','source','study','tissue','info'])
GeneSNP_Association = collections.namedtuple('GeneSNP_Association', ['gene', 'snp', 'score', 'rank', 'intermediary_scores', 'cisregulatory_evidence', 'regulatory_evidence'])
GeneCluster_Association = collections.namedtuple('GeneCluster_Association', ['gene', 'cluster', 'score', 'evidence','r2'])
