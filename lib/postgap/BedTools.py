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
import os.path
import pybedtools

def overlap_snps_to_bed(snps, bed):
	'''
		Find overlaps between SNP elements and annotated Bed file
		Args:
		* [ SNP ]
		* string (location of bed file)
		Returntype: pybedtools Interval iterator

	'''
	snps = list(snps)
	
	if (len(snps)==0):
		return []
	
	max_pos = max([x.pos for x in snps])
	min_pos = min([x.pos-1 for x in snps])
	chrom = snps[0].chrom
	SNP_bt = snps_to_bt(snps)
	Annotation_bt_indexed = bed_to_bt_indexed(bed)

	intersection = Annotation_bt_indexed.tabix_intervals('{}:{}-{}'.format(chrom,min_pos,max_pos)).intersect(SNP_bt, wa=True, wb=True)
	return intersection

def closest(snps, bed):
	SNP_string = "\n".join("\t".join((snp.chrom, str(snp.pos-1), str(snp.pos), snp.rsID)) for snp in sorted(snps, key=lambda X: (X.chrom, X.pos)))
	SNP_bt = pybedtools.BedTool(SNP_string, from_string=True)
	Annotation_bt_indexed = bed_to_bt_indexed(bed)
	return SNP_bt.closest(Annotation_bt_indexed, wa=True, wb=True)

def snps_to_bt(snps):
	SNP_string = "\n".join(( "\t".join((snp.chrom, str(snp.pos-1), str(snp.pos), snp.rsID)) for snp in snps ))
	return pybedtools.BedTool(SNP_string, from_string=True)

def bed_to_bt_indexed(bed):
	if not os.path.isfile(bed + '.gz') or not os.path.isfile(bed + '.gz.tbi'):
		Annotation_bt = pybedtools.BedTool(bed)
		return Annotation_bt.tabix()
	else:
		return pybedtools.BedTool(bed + '.gz')
